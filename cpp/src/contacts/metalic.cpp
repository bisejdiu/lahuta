#include "contacts/metals.hpp"
#include "lahuta.hpp"

namespace lahuta {

void find_metalic(const Luni *luni, GeometryOptions opts, Contacts &contacts) {

  double metal_distmax_ = 3.0;
  AtomDataVec metals = get_atom_data(luni, AtomType::IonicTypeMetal | AtomType::TransitionMetal);
  AtomDataVec metal_binders = get_atom_data(luni, AtomType::IonicTypePartner | AtomType::DativeBondPartner);

  auto m_grid = FastNS(metal_binders.positions(), metal_distmax_);
  auto m_nbrs = m_grid.search(metals.positions());

  for (const auto &[pair, dist] : m_nbrs) {
    auto [metal_index, metal_binding_index] = pair;
    const auto &metal = metals.data[metal_index];
    const auto &metal_binding = metal_binders.data[metal_binding_index];

    // NOTE: `is_metal_coordination` also checks for transition metal - transition metal
    // coordination, but our approach does not capture this interaction.
    // I think it is fine to ignore this case.
    // FIX: Given the type selection above, it's likely we do not need this check
    if (!is_metal_coordination(metal.type, metal_binding.type)
        && !is_metal_coordination(metal_binding.type, metal.type)) {
      std::cout << "Not a metal coordination bond" << std::endl;
      continue;
    }

    contacts.add(Contact(
        EntityID(metal.atom->getIdx()),
        EntityID(metal_binding.atom->getIdx()),
        dist,
        InteractionType::MetalCoordination));
  }
}

} // namespace lahuta
