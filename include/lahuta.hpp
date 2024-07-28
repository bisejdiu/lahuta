#include "bond_order.hpp"
#include "bonds.hpp"
#include "conv.hpp"
#include "nsgrid.hpp"
#include <optional>
#include <vector>

namespace Lahuta {

class ISource {
public:
  // FIX: RDKit provides ROMol and RWMol 
  virtual RDKit::RWMol &getMolecule() = 0;
  virtual const RDKit::RWMol &getMolecule() const = 0;
  virtual RDKit::Conformer &getConformer(int id = -1) = 0;
  virtual ~ISource() = default;

  virtual void process(const gemmi::Structure &st) = 0;
};

class GemmiSource : public ISource {
private:
  std::unique_ptr<RDKit::RWMol> mol = std::make_unique<RDKit::RWMol>();

public:
  explicit GemmiSource() = default;

  void process(const gemmi::Structure &st) override {
    RDKit::Conformer *conformer = new RDKit::Conformer();
    gemmiStructureToRDKit(*mol, st, *conformer, false);
    mol->addConformer(conformer, true);
  }

  RDKit::RWMol &getMolecule() override { return *mol; }
  const RDKit::RWMol &getMolecule() const override { return *mol; }

  RDKit::Conformer &getConformer(int id = -1) override {
    return mol->getConformer(id);
  }

};

class BondComputation {
public:
  static RDKit::RWMol computeProteinBonds(RDKit::RWMol &mol,
                                          const NSResults &neighbors) {
    std::vector<int> non_protein_indices;
    return lahutaBondAssignment(mol, neighbors, non_protein_indices);
  }

  static void perceiveBondOrders(RDKit::RWMol &mol) {
    mol.updatePropertyCache(false);
    PerceiveBondOrders(mol);
  }

  static void mergeBonds(RDKit::RWMol &targetMol, RDKit::RWMol &sourceMol,
                         const std::vector<int> &indexMap) {
    for (auto bondIt = sourceMol.beginBonds(); bondIt != sourceMol.endBonds();
         ++bondIt) {
      const RDKit::Bond *bond = *bondIt;
      int bIdx = indexMap[bond->getBeginAtomIdx()];
      int eIdx = indexMap[bond->getEndAtomIdx()];
      if (targetMol.getBondBetweenAtoms(bIdx, eIdx) == nullptr) {
        targetMol.addBond(bIdx, eIdx, bond->getBondType());
      }
    }
  }

  static void computeBonds(RDKit::RWMol &mol,
                           const NSResults &neighborResults) {

    // RDKit::RWMol resultMol = mol;
    // auto proteinMol = computeProteinBonds(resultMol, neighborResults);

    std::vector<int> non_protein_indices;
    auto newMol =
        lahutaBondAssignment(mol, neighborResults, non_protein_indices);
    newMol.updatePropertyCache(false);
    perceiveBondOrders(newMol);

    for (auto bondIt = newMol.beginBonds(); bondIt != newMol.endBonds();
         ++bondIt) {
      RDKit::Bond *bond = *bondIt;
      auto bAtomIdx = bond->getBeginAtomIdx();
      auto eAtomIdx = bond->getEndAtomIdx();
      int bIdx = non_protein_indices[bAtomIdx];
      int eIdx = non_protein_indices[eAtomIdx];
      if (mol.getBondBetweenAtoms(bIdx, eIdx) == nullptr) {
        mol.addBond(bIdx, eIdx, bond->getBondType());
      }
    }
  }
};

class Luni {
private:
  ISource &source;
  NSResults neighborResults;
  double _cutoff;

public:
  explicit Luni(ISource &source) : source(source), _cutoff(4.5) {
    initAndCompBonds();
  }

  void initAndCompBonds() {
    const auto &conf = source.getConformer();
    FastNS grid(conf.getPositions(), _cutoff);
    neighborResults = grid.selfSearch();
    BondComputation::computeBonds(source.getMolecule(), neighborResults);
  }

  RDKit::RWMol &getMolecule() { return source.getMolecule(); }
  const RDKit::RWMol &getMolecule() const { return source.getMolecule(); }
  const std::vector<RDGeom::Point3D> &getPositions() const {
    return source.getConformer().getPositions();
  }
  const NSResults &getNeighborResults() const { return neighborResults; }
  double getCutoff() const { return _cutoff; }

  void setNeighborResults(const NSResults &results) {
    neighborResults = results;
  }

  NSResults findNeighbors(std::optional<double> cutoff) const {
    auto value = cutoff.value_or(_cutoff);

    if (value == _cutoff) {
      return neighborResults;
    } else if (value < _cutoff) {
      return neighborResults.filterByDistance(value);
    }

    FastNS grid(getPositions(), value);
    return grid.selfSearch();
  }
};

} // namespace Lahuta
