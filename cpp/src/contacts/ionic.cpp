#include "contacts/ionic.hpp"
#include "contacts/search.hpp"
#include "lahuta.hpp"

namespace lahuta {

void find_ionic(const Luni &luni, GeometryOptions opts, Contacts &container) {

  const auto &conf = luni.get_molecule().getConformer();

  const FeatureVec positives = get_features(&luni, AtomType::POS_IONISABLE);
  const FeatureVec negatives = get_features(&luni, AtomType::NEG_IONISABLE);

  if (positives.get_data().empty() || negatives.get_data().empty()) {
    return;
  }

  double max_dist = 5.0;
  EntityNeighborSearch ens(conf);
  auto results = ens.search(positives, negatives, max_dist);

  for (const auto &[pair, dist] : results) {
    auto [positive_idx, negative_idx] = pair;
    auto positive = positives[positive_idx];
    auto negative = negatives[negative_idx];

    EntityID entity1 = make_entity_id(EntityType::Group, negative.get_id());
    EntityID entity2 = make_entity_id(EntityType::Group, positive.get_id());

    container.add(Contact(entity1, entity2, dist, InteractionType::Ionic));
  }
}

} // namespace lahuta
