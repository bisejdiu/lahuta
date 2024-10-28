#include "contacts/groups.hpp"
#include "contacts/aromaticity.hpp"
#include "contacts/charges.hpp"
#include "contacts/features.hpp"

namespace lahuta {

FeatureVec PositiveChargeGroup::identify(const RDKit::RWMol &mol, Residues &residues) const {
  return add_positive_charges(mol, residues);
};

FeatureVec NegativeChargeGroup::identify(const RDKit::RWMol &mol, Residues &residues) const {
  return add_negative_charges(mol, residues);
};

FeatureVec AromaticRingGroup::identify(const RDKit::RWMol &mol, Residues &residues) const {
  if (mol.getRingInfo()->isInitialized()) {
    mol.getRingInfo()->reset();
  }
  mol.getRingInfo()->initialize(RDKit::FIND_RING_TYPE_SYMM_SSSR);
  return add_aromatic_rings(mol, residues);
};

FeatureVec GroupTypeStrategy::identify(const RDKit::RWMol &mol) const {
  auto residues = Residues(mol);

  std::vector<Feature> group_features; // all features from all strategies
  for (const auto &strategy : strategies) {
    auto features = strategy->identify(mol, residues); // features from one strategy
    group_features.insert(group_features.end(), features.features.begin(), features.features.end());
  }
  // NOTE: `group_features` is not sorted as we get it from the strategies. If we do not sort `members`
  // somehow, the currently assign ids just lock-in the order of the features as it already is.
  assign_ids(group_features);
  compute_centers(group_features);

  return {group_features};
}

void GroupTypeStrategy::assign_ids(std::vector<Feature> &features) const {
  for (size_t i = 0; i < features.size(); ++i) {
    features[i].id = i;
  }
}

void GroupTypeStrategy::compute_centers(std::vector<Feature> &features) const {
  auto &conf = features.front().members.front()->getOwningMol().getConformer();
  for (auto &feature : features) {
    RDGeom::Point3D center_ = {0.0, 0.0, 0.0};
    for (const auto *atom : feature.members) {
      auto pos = conf.getAtomPos(atom->getIdx());
      center_ += pos;
    }
    center_ /= feature.members.size();
    feature.center[0] = center_.x;
    feature.center[1] = center_.y;
    feature.center[2] = center_.z;
  }
}

std::vector<const Feature *> get_features(const std::vector<Feature> &features, AtomType type) {
  std::vector<const Feature *> filtered_features;
  for (const auto &feature : features) {
    if (feature.type == type) {
      filtered_features.push_back(&feature);
    }
  }
  return filtered_features;
}

} // namespace lahuta
