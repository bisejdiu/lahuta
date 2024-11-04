#include "contacts/pistacking.hpp"
#include "contacts/geometry.hpp"
#include "contacts/search.hpp"
#include "lahuta.hpp"

namespace lahuta {

Contacts find_pistacking(const Luni &luni, PiStackingParams opts) {

  Contacts contacts(&luni);
  const auto rings = luni.get_rings();

  double dist_max = 6.0;
  EntityNeighborSearch ens(luni.get_conformer());
  auto nbrs = ens.search(rings, dist_max);

  for (const auto &[pair, dist] : nbrs) {
    auto [ring_index_a, ring_index_b] = pair;
    const auto &ring_a = rings[ring_index_a];
    const auto &ring_b = rings[ring_index_b];

    // e.g. trp
    if (is_same_residue(luni.get_molecule(), *ring_a.atoms.front(), *ring_b.atoms.front())) continue;

    auto dot_product = ring_a.norm.dotProduct(ring_b.norm);
    auto angle = std::acos(std::clamp(dot_product, -1.0, 1.0));
    if (angle > M_PI / 2) angle = M_PI - angle; // obtuse -> acute

    double offset_a = geometry::compute_in_plane_offset(ring_a.center, ring_b.center, ring_a.norm);
    double offset_b = geometry::compute_in_plane_offset(ring_b.center, ring_a.center, ring_b.norm);

    if (std::min(offset_a, offset_b) <= opts.offset_max) {
      EntityID entity1 = make_entity_id(EntityType::Ring, ring_a.get_id());
      EntityID entity2 = make_entity_id(EntityType::Ring, ring_b.get_id());

      (angle <= opts.angle_dev_max)
          ? contacts.add(Contact(entity1, entity2, dist, InteractionType::PiStackingP))
      : (std::abs(angle - M_PI / 2) <= opts.angle_dev_max)
          ? contacts.add(Contact(entity1, entity2, dist, InteractionType::PiStackingT))
          : void();
    }
  }

  return contacts;
}

} // namespace lahuta
