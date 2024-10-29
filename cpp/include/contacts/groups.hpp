#ifndef LAHUTA_GROUPS_HPP
#define LAHUTA_GROUPS_HPP

#include "GraphMol/RWMol.h"
#include "atom_types.hpp"
#include "contacts/aromaticity.hpp"
#include "contacts/charges.hpp"
#include "features.hpp"
#include "residues.hpp"

namespace lahuta {

class GroupTypeBase {
public:
  virtual ~GroupTypeBase() = default;
  virtual FeatureVec identify(const RDKit::RWMol &mol, const Residues &residues) const = 0;
};

class PositiveChargeGroup : public GroupTypeBase {
public:
  FeatureVec identify(const RDKit::RWMol &mol, const Residues &residues) const override {
    return add_positive_charges(mol, residues);
  }
};

class NegativeChargeGroup : public GroupTypeBase {
public:
  FeatureVec identify(const RDKit::RWMol &mol, const Residues &residues) const override {
    return add_negative_charges(mol, residues);
  };
};

class AromaticRingGroup : public GroupTypeBase {
public:
  FeatureVec identify(const RDKit::RWMol &mol, const Residues &residues) const override {
    return add_aromatic_rings(mol, residues);
  };
};

class GroupTypeStrategy {
private:
  std::vector<std::unique_ptr<GroupTypeBase>> strategies;

public:
  template <typename T> void add_strategy() {
    static_assert(std::is_base_of<GroupTypeBase, T>::value, "T must derive from GroupTypeStrategy");
    strategies.push_back(std::make_unique<T>());
  }

  FeatureVec identify(const RDKit::RWMol &mol, const Residues &residues) const;

private:
  void assign_ids(std::vector<Feature> &features) const;
  void compute_centers(std::vector<Feature> &features) const;
};

class GroupTypeFactory {
public:
  static GroupTypeStrategy create() {
    GroupTypeStrategy composite;
    composite.add_strategy<PositiveChargeGroup>();
    composite.add_strategy<NegativeChargeGroup>();
    composite.add_strategy<AromaticRingGroup>();
    return composite;
  }
};

class GroupTypeAnalysis {
public:
  static FeatureVec analyze(const RDKit::RWMol &mol, const Residues &residues) {
    auto strategy = GroupTypeFactory::create();
    return strategy.identify(mol, residues);
  }
};

std::vector<const Feature *> get_features(const std::vector<Feature> &features, AtomType type);

} // namespace lahuta

#endif // LAHUTA_GROUPS_HPP
