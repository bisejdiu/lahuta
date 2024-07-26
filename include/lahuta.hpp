#include "bonds.hpp"  // for perceiveBonds
#include "conv.hpp"   // for gemmiStructureToRDKit
#include "nsgrid.hpp" // for FastNS

#include "bond_order.hpp"

using namespace gemmi;
using namespace RDKit;

struct Lahuta {
private:
  RDKit::RWMol mol;
  std::vector<RDGeom::Point3D> atom_coords;
  NSResults results;

public:
  Lahuta(Structure &st) {
    RDKit::Conformer *conf = new RDKit::Conformer();
    RDKit::RWMol mol = gemmiStructureToRDKit(st, *conf, false);
    mol.addConformer(conf, true);

    double cutoff = 4.5;
    std::vector<RDGeom::Point3D> atom_coords = conf->getPositions();
    FastNS grid(atom_coords, cutoff);
    auto results = grid.selfSearch();

    std::vector<int> non_protein_indices;
    auto newMol = lahutaBondAssignment(mol, results, non_protein_indices);

    auto newMolConf = newMol.getConformer();
    newMol.updatePropertyCache(false);
    PerceiveBondOrders(newMol);

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

  RDKit::RWMol getMol() { return mol; }

  std::vector<RDGeom::Point3D> getAtomCoords() { return atom_coords; }

  NSResults getNSResult() { return results; }
};
