#include "contacts/ionic.hpp"
#include "contacts/search.hpp"
#include "lahuta.hpp"

namespace lahuta {

Contacts find_ionic(const Luni &luni, IonicParams opts) {

  Contacts contacts(&luni);
  const auto positives = GroupEntityCollection::filter(&luni, AtomType::POS_IONISABLE);
  const auto negatives = GroupEntityCollection::filter(&luni, AtomType::NEG_IONISABLE);

  EntityNeighborSearch ens(luni.get_conformer());
  auto results = ens.search(positives, negatives, opts.distance_max);

  for (const auto &[pair, dist] : results) {
    auto [positive_idx, negative_idx] = pair;
    auto positive = positives[positive_idx];
    auto negative = negatives[negative_idx];

    EntityID entity1 = make_entity_id(EntityType::Group, negative.get_id());
    EntityID entity2 = make_entity_id(EntityType::Group, positive.get_id());

    contacts.add(Contact(entity1, entity2, dist, InteractionType::Ionic));
  }

  return contacts;
}

} // namespace lahuta
