#include <gemmi/mmread_gz.hpp> // for read_structure_gz

#include "bonds.hpp"  // for perceiveBonds
#include "conv.hpp"   // for gemmiStructureToRDKit
#include "nsgrid.hpp" // for FastNS

#include <chrono>

#include "bond_order.hpp"
#include "bond_utils.hpp" // for CleanUpMolecule

#include "bond_table/bonds.hpp"
#include "bond_table/table.hpp"

#include "kekulize.h"

#include <iostream>

#define TO_MS(d) std::chrono::duration_cast<std::chrono::milliseconds>(d)
#define T() std::chrono::high_resolution_clock::now()

using namespace gemmi;
using namespace RDKit;

int main(int argc, char const *argv[]) {
  auto startTotTime = std::chrono::high_resolution_clock::now();
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <cif file>" << std::endl;
    return 1;
  }
  std::string file_name = argv[1];

  Structure st = read_structure_gz(file_name);
  std::cout << "Done: Loading structure" << std::endl;

  RDKit::Conformer *conf = new RDKit::Conformer();
  RDKit::RWMol mol = gemmiStructureToRDKit(st, *conf, false);
  mol.addConformer(conf, true);

  double cutoff = 4.5; // 01;
  std::vector<RDGeom::Point3D> atom_coords = conf->getPositions();
  FastNS grid(atom_coords, cutoff);
  auto results = grid.selfSearch();
  std::cout << "Done: Computing Neighbors" << results.getNeighborPairsSize()
            << std::endl;

  // if (!flag_b) {
  //   perceiveBonds(mol, results, 0.45);
  //
  //   mol.updatePropertyCache(false);
  //   CleanUpMolecule(mol, *conf);
  //
  //   PerceiveBondOrders(mol);
  // }

  auto start = T();
  std::vector<int> non_protein_indices;
  auto newMol = lahutaBondAssignment(mol, results, non_protein_indices);
  auto end = T();
  std::cout << "T1" << TO_MS(end - start).count() << "ms" << std::endl;

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

  // for (gemmi::Connection &conn : st.connections) {
  //   // FIX: Need to iterate over all models
  //   gemmi::Atom *a1 = st.first_model().find_cra(conn.partner1).atom;
  //   gemmi::Atom *a2 = st.first_model().find_cra(conn.partner2).atom;
  //
  //   if (mol.getBondBetweenAtoms(a1->serial - 1, a2->serial - 1) == nullptr) {
  //     // compute distance between the two atoms:
  //     double dist =
  //         (atom_coords[a1->serial - 1] - atom_coords[a2->serial -
  //         1]).length();
  //     std::cout << "conn pair: " << a1->serial << " " << a2->serial << " "
  //               << dist << std::endl;
  //     // mol.addBond((unsigned int)a1->serial - 1, (unsigned int)a2->serial -
  //     1,
  //     //             RDKit::Bond::BondType::SINGLE);
  //     // continue;
  //   }
  // }

  // std::string bondOrders = "";
  // for (auto bondIt = mol.beginBonds(); bondIt != mol.endBonds(); ++bondIt) {
  //   RDKit::Bond *bond = *bondIt;
  //   auto res1 = mol.getAtomWithIdx(bond->getBeginAtomIdx())->getMonomerInfo();
  //   auto res2 = mol.getAtomWithIdx(bond->getEndAtomIdx())->getMonomerInfo();
  //
  //   AtomPDBResidueInfo *residue1 = dynamic_cast<AtomPDBResidueInfo *>(res1);
  //   AtomPDBResidueInfo *residue2 = dynamic_cast<AtomPDBResidueInfo *>(res2);
  //
  //   bondOrders += std::to_string(bond->getBeginAtomIdx()) + " " +
  //                 // res1->getName() + " " + residue1->getResidueName() + " " +
  //                 // residue1->getAltLoc() + " " +
  //                 std::to_string(bond->getEndAtomIdx()) + " " +
  //                 // res2->getName() + " " + residue2->getResidueName() + " " +
  //                 // residue2->getAltLoc() + " " +
  //                 std::to_string(bond->getBondType()) + "\n";
  // }
  // std::cout << "FINAL RESULT: " << std::endl;
  // std::cout << bondOrders << std::endl;

  return 0;
}
