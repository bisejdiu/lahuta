#include "contacts/groups.hpp"

namespace lahuta {

GroupEntityCollection GroupTypeStrategy::identify(const RDKit::RWMol &mol, const Residues &residues) const {

  std::vector<GroupEntity> group_features; // all features from all strategies
  for (const auto &strategy : strategies) {
    auto features = strategy->identify(mol, residues); // features from one strategy
    group_features.insert(group_features.end(), features.get_data().begin(), features.get_data().end());
  }
  // NOTE: `group_features` is not sorted as we get it from the strategies.
  // If we do not sort `members` somehow, the currently assign ids just lock-in
  // the order of the features as it already is.
  assign_ids(group_features);
  compute_centers(group_features);

  return {group_features};
}

void GroupTypeStrategy::assign_ids(std::vector<GroupEntity> &features) const {
  for (size_t i = 0; i < features.size(); ++i) {
    features[i].set_id(i);
  }
}

void GroupTypeStrategy::compute_centers(std::vector<GroupEntity> &features) const {
  if (features.empty() || features.front().atoms.empty()) return;

  auto &conf = features.front().atoms.front()->getOwningMol().getConformer();
  for (auto &feature : features) {
    RDGeom::Point3D center_ = {0.0, 0.0, 0.0};
    for (const auto *atom : feature.atoms) {
      auto pos = conf.getAtomPos(atom->getIdx());
      center_ += pos;
    }
    center_ /= feature.atoms.size();
    feature.center[0] = center_.x;
    feature.center[1] = center_.y;
    feature.center[2] = center_.z;
  }
}

} // namespace lahuta
