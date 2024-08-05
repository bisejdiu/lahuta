#include <GraphMol/RDKitBase.h>

#include "bond_order.hpp"
#include "bonds.hpp"
#include "convert.hpp"
#include "nsgrid.hpp"
#include <optional>
#include <vector>
#include "ob/clean_mol.hpp"
#include <gemmi/mmread_gz.hpp> // for read_structure_gz
//
#include <chrono>

#define TO_MS(d) std::chrono::duration_cast<std::chrono::milliseconds>(d)
#define T() std::chrono::high_resolution_clock::now()

namespace Lahuta {

class ISource {
public:
  // FIX: RDKit provides ROMol and RWMol
  virtual RDKit::RWMol &getMolecule() = 0;
  virtual const RDKit::RWMol &getMolecule() const = 0;
  virtual RDKit::Conformer &getConformer(int id = -1) = 0;
  virtual ~ISource() = default;

  // virtual void process(const gemmi::Structure &st) = 0;
  virtual void process(std::string file_name) = 0;
};

class GemmiSource : public ISource {
private:
  std::shared_ptr<RDKit::RWMol> mol = std::make_shared<RDKit::RWMol>();

public:
  explicit GemmiSource() = default;

  // copy constructor
  GemmiSource(const GemmiSource &source) {
    mol = std::make_unique<RDKit::RWMol>(*source.mol);
  }

  // void process(const gemmi::Structure &st) override {
  void process(std::string file_name) override {
    RDKit::Conformer *conformer = new RDKit::Conformer();
    Structure st = read_structure_gz(file_name);
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
    std::vector<int> non_predef_atom_indices;
    return assign_bonds(mol, neighbors, non_predef_atom_indices);
  }

  static void perceiveBondOrders(RDKit::RWMol &mol) {
    mol.updatePropertyCache(false);
    perceive_bond_orders_obabel(mol);
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
    non_protein_indices.reserve(mol.getNumAtoms());
    auto newMol =
        assign_bonds(mol, neighborResults, non_protein_indices);
    // compute time to clean molecule in ms:
    auto start = T();
    clean_bonds(newMol, newMol.getConformer());
    auto cleanTime = TO_MS(T() - start).count();
    std::cout << "Time to clean molecule: " << cleanTime << "ms" << std::endl;


    
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
  std::unique_ptr<ISource> source;
  NSResults neighborResults;
  float _cutoff;
  FastNS grid;

public:
  // explicit Luni(GemmiSource &source) : source(&source), _cutoff(4.5) { 
  //   initAndCompBonds();
  // }
  explicit Luni(std::string file_name) : _cutoff(4.5) { 
    source = std::make_unique<GemmiSource>();
    source->process(file_name);
    initAndCompBonds();
  }

  void initAndCompBonds() {
    const auto &conf = source->getConformer();
    grid = FastNS(conf.getPositions(), _cutoff);
    neighborResults = grid.selfSearch();
    BondComputation::computeBonds(source->getMolecule(), neighborResults);
  }

  RDKit::RWMol &getMolecule() { return source->getMolecule(); }
  const RDKit::RWMol &getMolecule() const { return source->getMolecule(); }
  const std::vector<RDGeom::Point3D> &getPositions(int confId = -1) const {
    return source->getConformer(confId).getPositions();
  }
  const NSResults &getNeighborResults() const { return neighborResults; }
  double getCutoff() const { return _cutoff; }

  void setNeighborResults(const NSResults &results) {
    neighborResults = results;
  }

  NSResults findNeighbors(std::optional<double> cutoff) {
    auto value = cutoff.value_or(_cutoff);

    if (value == _cutoff) {
      return neighborResults;
    } else if (value < _cutoff) {
      return neighborResults.filterByDistance(value);
    }

    grid.updateCutoff(value);
    return grid.selfSearch();
  }
};

} // namespace Lahuta
