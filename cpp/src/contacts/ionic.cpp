#include "contacts/ionic.hpp"
#include "contacts/search.hpp"
#include "lahuta.hpp"
#include "common.hpp"

namespace lahuta {

Contacts find_ionic(const Luni &luni, IonicParams opts) {

  Contacts contacts(&luni);
  const auto positives = GroupEntityCollection::filter(&luni, AtomType::POS_IONISABLE);
  const auto negatives = GroupEntityCollection::filter(&luni, AtomType::NEG_IONISABLE);

  EntityNeighborSearch ens(luni.get_conformer());
  auto results = ens.search(positives, negatives, opts.distance_max);

  std::unordered_set<std::pair<int, int>, common::PairHash> seen;

  for (const auto &[pair, dist] : results) {
    auto [positive_idx, negative_idx] = pair;
    auto positive = positives[positive_idx];
    auto negative = negatives[negative_idx];

    if (positive.get_id() == negative.get_id() || dist < 2.0) continue;
    if (is_duplicate({negative.get_id(), positive.get_id()}, seen)) continue;

    EntityID entity1 = make_entity_id(EntityType::Group, negative.get_id());
    EntityID entity2 = make_entity_id(EntityType::Group, positive.get_id());

    contacts.add(Contact(entity1, entity2, dist, InteractionType::Ionic));
  }

  return contacts;
}

} // namespace lahuta
