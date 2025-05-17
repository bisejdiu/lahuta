#include "contacts/metalic.hpp"
#include "contacts/search.hpp"
#include "lahuta.hpp"

namespace lahuta {

bool is_metalic(AtomType at1, AtomType at2) {
  using AtomTypeFlags::has;
  if (has(at1, AtomType::TransitionMetal)) return has(at2, AtomType::DativeBondPartner);
  if (has(at1, AtomType::IonicTypeMetal))  return has(at2, AtomType::IonicTypePartner);
  return false;
}

Contacts find_metalic(const Luni &luni, std::optional<MetalicParams> params) {
  using AEC = AtomEntityCollection;

  Contacts contacts(&luni);
  const auto &mol = luni.get_molecule();
  MetalicParams opts = params.value_or(MetalicParams{});

  AtomEntityCollection metals, metal_binders;
  metals = AEC::filter(&luni, AtomType::IonicTypeMetal | AtomType::TransitionMetal);
  metal_binders = AEC::filter(&luni, AtomType::IonicTypePartner | AtomType::DativeBondPartner);

  auto m_nbrs = EntityNeighborSearch::search(metals, metal_binders, opts.distance_max);

  for (const auto &[pair, dist] : m_nbrs) {
    auto [metal_index, metal_binding_index] = pair;
    const auto &metal = metals.get_data()[metal_index];
    const auto &metal_binding = metal_binders.get_data()[metal_binding_index];

    if (!is_metalic(metal.type, metal_binding.type) && !is_metalic(metal_binding.type, metal.type)) continue;

    auto bond = mol.getBondBetweenAtoms(metal.atom->getIdx(), metal_binding.atom->getIdx());
    if (bond) continue;

    contacts.add(Contact(
        static_cast<EntityID>(metal.atom->getIdx()),
        static_cast<EntityID>(metal_binding.atom->getIdx()),
        dist,
        InteractionType::MetalCoordination));
  }

  return contacts;
}

} // namespace lahuta
