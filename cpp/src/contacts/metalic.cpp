#include "contacts/search.hpp"
#include "lahuta.hpp"

namespace lahuta {

bool is_metalic(AtomType at1, AtomType at2) {
  using AtomTypeFlags::has;
  if (has(at1, AtomType::TransitionMetal)) return has(at2, AtomType::DativeBondPartner);
  if (has(at1, AtomType::IonicTypeMetal)) return has(at2, AtomType::IonicTypePartner);
  return false;
}

void find_metalic(const Luni *luni, GeometryOptions opts, Contacts &contacts) {

  double metal_distmax_ = 3.0;
  AtomDataVec metals = get_atom_data(luni, AtomType::IonicTypeMetal | AtomType::TransitionMetal);
  AtomDataVec metal_binders = get_atom_data(luni, AtomType::IonicTypePartner | AtomType::DativeBondPartner);

  EntityNeighborSearch ens(luni->get_molecule().getConformer());
  auto m_nbrs = ens.search(metals, metal_binders, metal_distmax_);

  for (const auto &[pair, dist] : m_nbrs) {
    auto [metal_index, metal_binding_index] = pair;
    const auto &metal = metals.get_data()[metal_index];
    const auto &metal_binding = metal_binders.get_data()[metal_binding_index];

    if (!is_metalic(metal.type, metal_binding.type) && !is_metalic(metal_binding.type, metal.type)) {
      continue;
    }

    contacts.add(Contact(
        static_cast<EntityID>(metal.atom->getIdx()),
        static_cast<EntityID>(metal_binding.atom->getIdx()),
        dist,
        InteractionType::MetalCoordination));
  }
}

} // namespace lahuta
