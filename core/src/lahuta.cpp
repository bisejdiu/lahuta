#include <stdexcept>
#include <string>

#include "GraphMol/PeriodicTable.h"
#include "gemmi/gz.hpp"
#include "lahuta.hpp"
#include "mmap/MemoryMapped.h"
#include "models/factory.hpp"
#include "models/parser.hpp"
#include "models/topology.hpp"
#include "nsgrid.hpp"
#include "read_file.hpp"
#include "selections/mol_filters.hpp"

// clang-format off
namespace lahuta {

Luni::Luni(std::string file_name) : file_name_(file_name) {
  Logger::get_logger()->debug("Processing file: {}", file_name_);
  mol = read_and_make_molecule(gemmi::MaybeGzipped(file_name_));

  /*auto st = gemmi::read_structure_gz(file_name_);*/
  /*mol = std::make_shared<RDKit::RWMol>();*/
  /*RDKit::Conformer *conformer = new RDKit::Conformer();*/
  /*create_RDKit_repr(*mol, st, *conformer);*/
  /*mol->updatePropertyCache(false);*/
  /*mol->addConformer(conformer, true);*/
}

Luni Luni::create(const gemmi::Structure &st) {
  auto mol = create_RDKit(st);
  return Luni(mol);
}

Luni Luni::create(const IR &ir) {
  auto mol = std::make_shared<RDKit::RWMol>();
  IR_to_RWMol(*mol, ir);
  return Luni(mol);
}

Luni::Luni(std::string file_name, ModelFileTag) : file_name_(file_name) {
  // Ensure model pools are initialized for single threaded use by default.
  // The pipeline path will reinitialize with the correct thread count.
  try {
    // Pool initialization is idempotent
    InfoPoolFactory::initialize(1);
    BondPoolFactory::initialize(1);
    AtomPoolFactory::initialize(1);
  } catch (...) {}

  model_origin_ = true;
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
    if (!build_model_topology(mol, result, ModelTopologyMethod::CSR)) {
      throw std::runtime_error("Failed to build model topology from model file");
    }
  } catch (const std::exception &e) {
    Logger::get_logger()->critical("Exception processing file {}: {}", file_name_, e.what());
  } catch (...) {
    Logger::get_logger()->critical("Unknown exception processing file: {}", file_name_);
  }
}

bool Luni::build_topology(std::optional<TopologyBuildingOptions> tops) {
  try {
    if (topology_built_) {
      Logger::get_logger()->debug("Topology already built, skipping rebuild to prevent molecule state corruption");
      return true;
    }

    ensure_topology_initialized();

    if (tops) {
      topology_built_ = topology->build(*tops);
    } else {
      TopologyBuildingOptions opts;
      if (model_origin_) opts.mode = TopologyBuildMode::Model;
      topology_built_ = topology->build(opts);
    }

    return topology_built_;
  } catch (const std::exception &e) {
    Logger::get_logger()->error("Error building topology: {}", e.what());
    return false;
  }
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

} // namespace lahuta
