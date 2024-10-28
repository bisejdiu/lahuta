#include "contacts/features.hpp"
#include "lahuta.hpp"

namespace lahuta {

/*std::vector<const Feature *> get_features(const std::vector<Feature> &features, AtomType type) {*/
/*  std::vector<const Feature *> filtered_features;*/
/*  for (const auto &feature : features) {*/
/*    if (feature.type == type) {*/
/*      filtered_features.push_back(&feature);*/
/*    }*/
/*  }*/
/*  return filtered_features;*/
/*}*/

Feature create_feature(AtomType type, FeatureGroup group, const std::vector<const RDKit::Atom *> &members) {
  Feature feature;
  feature.type = type;
  feature.group = group;
  feature.members = members;
  return feature;
}


FeatureVec get_features(const Luni *luni, AtomType type, FeatureTypeCheckFunc check_func) {
  auto &features = luni->get_features();
  FeatureVec feature_vec;
  for (const auto &feature : features.features) {
    if (check_func(feature.type, type)) {
      // confirm members is defined:
      if (feature.members.empty()) {
        throw std::runtime_error("Feature members are empty");
      }
      if (feature.members.front() == nullptr) {
        throw std::runtime_error("Feature members are not defined");
      }
      feature_vec.features.emplace_back(
          Feature{feature.type, feature.group, feature.members, feature.center});
      feature_vec.features.back().id = feature.id;
    }
  }

  return feature_vec;
}

} // namespace lahuta
