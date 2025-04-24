#include "lahuta.hpp"
#include "GraphMol/PeriodicTable.h"
#include "gemmi/gz.hpp"
#include "mmap/MemoryMapped.h"
#include "models/parser.hpp"
#include "models/topology.hpp"
#include "nsgrid.hpp"
#include <string>

namespace lahuta {

Luni::Luni(std::string file_name) : file_name_(file_name) {
  Logger::get_logger()->info("Processing file: {}", file_name_);
  read_structure();
}

Luni Luni::create(const gemmi::Structure &st) {
  auto mol = std::make_shared<RDKit::RWMol>();
  RDKit::Conformer *conformer = new RDKit::Conformer();
  create_RDKit_repr(*mol, st, *conformer, false);
  mol->updatePropertyCache(false);
  mol->addConformer(conformer, true);
  return Luni(mol);
}

Luni Luni::create(const IR &ir) {
  auto mol = std::make_shared<RDKit::RWMol>();
  IR_to_RWMol(*mol, ir);
  return Luni(mol);
}

Luni::Luni(std::string file_name, bool test) : file_name_(file_name) {

  ModelParserResult result;

  try {
    gemmi::MaybeGzipped input_file(file_name_);

    if (input_file.is_compressed()) {
      gemmi::CharArray buffer = input_file.uncompress_into_buffer();
      result = parse_model(buffer.data(), buffer.size());

    } else {
      MemoryMapped mm(file_name_.c_str());

      if (!mm.isValid()) {
        Logger::get_logger()->critical("Error opening file: {}", file_name_);
        return;
      }

      const char *data = reinterpret_cast<const char *>(mm.getData());
      size_t size = static_cast<size_t>(mm.size());

      result = parse_model(data, size);
    }
    build_model_topology(mol, result, ModelTopologyMethod::CSR);
  } catch (const std::exception &e) {
    Logger::get_logger()->critical("Exception processing file {}: {}", file_name_, e.what());
  } catch (...) {
    Logger::get_logger()->critical("Unknown exception processing file: {}", file_name_);
  }
}

bool Luni::build_topology(std::optional<TopologyBuildingOptions> tops) {
  try {
    ensure_topology_initialized();
    
    // If we're in filtered state, disable bond computation and dependent computations
    if (is_in_filtered_state) {
      topology->enable_only(TopologyComputation::None);
    }
    
    // Build the topology with the provided or default options
    if (tops) {
      topology->build(*tops);
    } else {
      // Use default options
      auto default_options = TopologyBuildingOptions{};
      topology->build(default_options);
    }
    
    // Mark topology as fully built
    topology_built_ = true;
    
    return true;
  } catch (const std::exception &e) {
    // Log the error but don't re-throw - just return false to indicate failure
    Logger::get_logger()->error("Error building topology: {}", e.what());
    
    // Still mark the topology as "built" but empty - this helps prevent accessing
    // non-existent data later
    topology_built_ = true;
    
    return false;
  }
}

void Luni::read_structure() {
  auto start1 = std::chrono::high_resolution_clock::now();
  auto st = gemmi::read_structure_gz(file_name_);
  auto end1 = std::chrono::high_resolution_clock::now();
  auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>(end1 - start1);
  Logger::get_logger()->info("gemmi Read file in {} us", duration1.count());

  auto start2 = std::chrono::high_resolution_clock::now();
  RDKit::Conformer *conformer = new RDKit::Conformer();
  create_RDKit_repr(*mol, st, *conformer, false);
  mol->updatePropertyCache(false);
  mol->addConformer(conformer, true);
  auto end2 = std::chrono::high_resolution_clock::now();
  auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(end2 - start2);
  Logger::get_logger()->info("RDKit Created mol in {} us", duration2.count());
}

const Topology *Luni::get_topology_ptr() const {
  if (!topology) {
    Logger::get_logger()->error("Cannot return topology, because no topology has been built.");
    throw std::runtime_error("Topology not built");
  }
  return &topology.value();
}

size_t Luni::total_size() const {
  size_t size = sizeof(Luni);

  size += mol->getNumAtoms() * sizeof(RDKit::Atom) + mol->getNumAtoms() * sizeof(RDKit::Atom *);
  size += mol->getNumBonds() * sizeof(RDKit::Bond) + mol->getNumBonds() * sizeof(RDKit::Bond *);
  size += mol->getNumConformers() * sizeof(RDKit::Conformer);
  size += mol->getNumConformers() * mol->getNumAtoms() * sizeof(RDGeom::Point3D);

  if (topology) {
    size += topology->total_size();
    size += mol->getRingInfo()->getTotalMemory();
  }

  size += entities.size() * sizeof(EntityType);
  for (const auto &[type, ids] : entities) {
    size += ids.size() * sizeof(EntityID);
  }

  size += sizeof(char) * file_name_.size();
  size += sizeof(int)  * filtered_indices.size();
  size += sizeof(bool);

  return size;
}

const std::vector<std::string> Luni::symbols() const {
  return atom_attrs<std::string>([](const RDKit::Atom *atom) { return atom->getSymbol(); });
}

const std::vector<std::string> Luni::names() const {
  return atom_attrs<std::string>([](const RDKit::Atom *atom) -> const std::string & {
    auto *info = static_cast<const RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
    return info->getName();
  });
}

const std::vector<int> Luni::indices() const {
  return atom_attrs<int>([](const RDKit::Atom *atom) { return atom->getIdx(); });
}

const std::vector<int> Luni::atomic_numbers() const {
  return atom_attrs<int>([](const RDKit::Atom *atom) { return atom->getAtomicNum(); });
}

const std::vector<std::string> Luni::elements() const {
  const RDKit::PeriodicTable *tbl = RDKit::PeriodicTable::getTable();
  return atom_attrs<std::string>(
      [&tbl](const RDKit::Atom *atom) { return tbl->getElementSymbol(atom->getAtomicNum()); });
}

const std::vector<std::string> Luni::resnames() const {
  return atom_attrs<std::string>([](const RDKit::Atom *atom) -> const std::string & {
    auto *info = static_cast<const RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
    return info->getResidueName();
  });
}

const std::vector<int> Luni::resids() const {
  return atom_attrs<int>([](const RDKit::Atom *atom) -> int {
    auto *info = static_cast<const RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
    return info->getResidueNumber();
  });
}

const std::vector<int> Luni::resindices() const {
  return atom_attrs<int>([](const RDKit::Atom *atom) -> int {
    auto *info = static_cast<const RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
    return info->getSegmentNumber();
  });
}

const std::vector<std::string> Luni::chainlabels() const {
  return atom_attrs<std::string>([](const RDKit::Atom *atom) -> const std::string & {
    auto *info = static_cast<const RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
    return info->getChainId();
  });
}

template <typename T> std::vector<T> Luni::atom_attrs(std::function<T(const RDKit::Atom *)> func) const {
  std::vector<T> attrs;
  attrs.reserve(mol->getNumAtoms());
  for (const auto atom : mol->atoms()) {
    attrs.push_back(func(atom));
  }
  return attrs;
}

template <typename T>
std::vector<std::reference_wrapper<const T>>
Luni::atom_attrs_ref(std::function<const T &(const RDKit::Atom *)> func) const {
  std::vector<std::reference_wrapper<const T>> attributes;
  attributes.reserve(mol->getNumAtoms());
  for (const auto atom : mol->atoms()) {
    attributes.push_back(std::cref(func(atom)));
  }
  return attributes;
}

auto Luni::match_smarts_string(std::string sm, std::string atype, bool log_values) const {

  if (!mol->getRingInfo()->isInitialized()) {
    RDKit::MolOps::symmetrizeSSSR(*mol);
  }
  std::vector<RDKit::MatchVectType> match_list;
  auto sm_mol = RDKit::SmartsToMol(sm);
  mol->updatePropertyCache(false);
  RDKit::SubstructMatch(*mol, *sm_mol, match_list);
  return match_list;
};

Luni Luni::filter(std::vector<int> &atom_indices) const {

  auto new_mol = filter_with_bonds(*mol, atom_indices);
  new_mol.updatePropertyCache(false);

  auto new_luni = Luni::create(std::make_shared<RDKit::RWMol>(new_mol));
  new_luni.is_in_filtered_state = true;
  return new_luni;
}

//
// The way Entities are generated, outside of atom entities, we don't have a guarantee that the index matches
// the entity id. If we decide to use these implementations, we need to use `get_*_entities` methods to first
// populate the entities. Only then can we use the `get_entity` method.  - Besian, March, 2025
//
template <> const RDKit::Atom &Luni::get_entity<RDKit::Atom>(EntityID id) const {
  auto index = get_entity_index(id);
  auto r = mol->getAtomWithIdx(index);
  return *r;
}

template <> const AtomEntity &Luni::get_entity<AtomEntity>(EntityID id) const {
  auto index = get_entity_index(id);
  return get_atom_types().get_data()[index];
}

template <> const RingEntity &Luni::get_entity<RingEntity>(EntityID id) const {
  auto index = get_entity_index(id);
  return get_rings().get_data()[index];
}

template <> const GroupEntity &Luni::get_entity<GroupEntity>(EntityID id) const {
  auto index = get_entity_index(id);
  // FIX: check if the index is valid
  return get_features()[index];
}

const std::vector<EntityID> &Luni::get_or_create_atom_entities() {
  if (entities.find(EntityType::Atom) == entities.end()) {
    std::vector<EntityID> atom_entities;
    auto mol = get_molecule();

    atom_entities.reserve(mol.getNumAtoms());
    for (const auto &atom : mol.atoms()) {
      atom_entities.push_back(make_entity_id(EntityType::Atom, atom->getIdx()));
    }
    entities[EntityType::Atom] = std::move(atom_entities);
  }

  return entities[EntityType::Atom];
}

const std::vector<EntityID> &Luni::get_or_create_ring_entities() {
  if (entities.find(EntityType::Ring) == entities.end()) {

    std::vector<EntityID> ring_entities;
    std::vector<RingEntity> rings = get_rings().get_data();

    ring_entities.reserve(rings.size());
    for (std::size_t i = 0; i < rings.size(); ++i) {
      ring_entities.push_back(make_entity_id(EntityType::Ring, i));
    }
    entities[EntityType::Ring] = std::move(ring_entities);
  }

  return entities[EntityType::Ring];
}

const std::vector<EntityID> &Luni::get_or_create_group_entities() {
  if (entities.find(EntityType::Group) == entities.end()) {
    std::vector<EntityID> group_entities;
    auto groups = get_features();
    group_entities.reserve(groups.get_data().size());
    for (std::size_t i = 0; i < groups.get_data().size(); ++i) {
      group_entities.push_back(make_entity_id(EntityType::Group, i));
    }
    entities[EntityType::Group] = std::move(group_entities);
  }

  return entities[EntityType::Group];
}

} // namespace lahuta
