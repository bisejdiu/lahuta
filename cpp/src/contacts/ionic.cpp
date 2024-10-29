#include "contacts/ionic.hpp"
#include "lahuta.hpp"

namespace lahuta {

void find_ionic(const Luni &luni, GeometryOptions opts, Contacts &container) {

  const auto &conf = luni.get_molecule().getConformer();

  const FeatureVec features_a = get_features(&luni, AtomType::POS_IONISABLE);
  const FeatureVec features_b = get_features(&luni, AtomType::NEG_IONISABLE);

  if (features_a.features.empty() || features_b.features.empty()) {
    return;
  }

  double max_dist = 5.0;
  FastNS grid(features_a.positions(), max_dist);
  auto atom_pairs = grid.search(features_b.positions());

  std::set<std::pair<size_t, size_t>> contacts;
  for (const auto &[pair, dist] : atom_pairs) {
    auto [feature_b_ix, feature_a_ix] = pair;
    auto feature_a = features_a[feature_a_ix];
    auto feature_b = features_b[feature_b_ix];

    auto feature_pair = std::pair{std::minmax(feature_b.get_id(), feature_a.get_id())};
    if (contacts.find(feature_pair) == contacts.end()) {
      contacts.insert(feature_pair);

      EntityID entity1 = make_entity_id(lahuta::EntityType::Group, feature_a.get_id());
      EntityID entity2 = make_entity_id(lahuta::EntityType::Group, feature_b.get_id());

      container.add(Contact(entity1, entity2, dist, InteractionType::Ionic));
    }
  }
}

} // namespace lahuta
