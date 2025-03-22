#ifndef LAHUTA_GROUPS_HPP
#define LAHUTA_GROUPS_HPP

#include "GraphMol/RWMol.h"
#include "contacts/aromaticity.hpp"
#include "contacts/charges.hpp"
#include "residues.hpp"

namespace lahuta {

class GroupTypeBase {
public:
  virtual ~GroupTypeBase() = default;
  virtual std::string name() const = 0;
  virtual GroupEntityCollection identify(const RDKit::RWMol &mol, const Residues &residues) const = 0;
};

class PositiveChargeGroup : public GroupTypeBase {
public:
  std::string name() const override { return "PositiveChargeGroup"; }
  GroupEntityCollection identify(const RDKit::RWMol &mol, const Residues &residues) const override {
    return add_positive_charges(mol, residues);
  }
};

class NegativeChargeGroup : public GroupTypeBase {
public:
  std::string name() const override { return "NegativeChargeGroup"; }
  GroupEntityCollection identify(const RDKit::RWMol &mol, const Residues &residues) const override {
    return add_negative_charges(mol, residues);
  };
};

class AromaticRingGroup : public GroupTypeBase {
public:
  std::string name() const override { return "AromaticRingGroup"; }
  GroupEntityCollection identify(const RDKit::RWMol &mol, const Residues &residues) const override {
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

  GroupEntityCollection identify(const RDKit::RWMol &mol, const Residues &residues) const;

private:
  void assign_ids(std::vector<GroupEntity> &features) const;
  void compute_centers(std::vector<GroupEntity> &features) const;
};

class GroupTypeFactory {
public:
  static GroupTypeStrategy create() {
    GroupTypeStrategy gt;
    // FIX: rename to add_method or similar, not add_strategy
    gt.add_strategy<PositiveChargeGroup>();
    gt.add_strategy<NegativeChargeGroup>();
    gt.add_strategy<AromaticRingGroup>();
    return gt;
  }
};

class GroupTypeAnalysis {
public:
  static GroupEntityCollection analyze(const RDKit::RWMol &mol, const Residues &residues) {
    auto strategy = GroupTypeFactory::create();
    return strategy.identify(mol, residues);
  }
};

} // namespace lahuta

#endif // LAHUTA_GROUPS_HPP
