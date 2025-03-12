#include "lahuta.hpp"
#include "nsgrid.hpp"

namespace lahuta {

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
// The way Entities are generated, outside of atom entities, we don't have a guarantee that the index matches the entity id.
// If we decide to use these implementations, we need to use `get_*_entities` methods to first populate the entities.
// Only then can we use the `get_entity` method.  - Besian, March, 2025
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
