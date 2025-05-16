#include "contacts/cationpi.hpp"
#include "contacts/search.hpp"
#include "contacts/utils.hpp"
#include "lahuta.hpp"
#include "neighbors.hpp"
#include "contacts/geometry.hpp"

namespace lahuta {

Contacts find_cationpi(const Luni &luni, CationPiParams opts) {

  Contacts contacts(&luni);

  const auto rings = luni.get_rings();
  const auto features = GroupEntityCollection::filter(&luni, AtomType::POS_IONISABLE);

  auto nbrs = EntityNeighborSearch::search(features, rings, opts.distance_max);

  for (const auto &[pair, dist] : nbrs) {
    auto [feature_index, ring_index] = pair;
    const auto &ring = rings.data[ring_index];
    const auto &feature = features[feature_index];

    auto first_ring_atom = ring.atoms.front();

    if (is_same_residue(luni.get_molecule(), *first_ring_atom, *feature.atoms.front())) continue;

    auto offset = geometry::compute_in_plane_offset(feature.center, ring.center, ring.normal);

    if (offset <= opts.offset_max) {
      contacts.add(Contact(
          make_entity_id(EntityType::Group, feature.get_id()),
          make_entity_id(EntityType::Ring, ring.get_id()),
          dist,
          InteractionType::CationPi));
    }
  }

  return contacts;
}

} // namespace lahuta
