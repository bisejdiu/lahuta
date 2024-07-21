#include <gemmi/mmread_gz.hpp> // for read_structure_gz

#include "bonds.hpp"  // for perceiveBonds
#include "conv.hpp"   // for gemmiStructureToRDKit
#include "nsgrid.hpp" // for FastNS

#include <chrono>

#include "bond_order.hpp"
#include "bond_utils.hpp" // for CleanUpMolecule

#include "bond_table/bonds.hpp"
#include "bond_table/table.hpp"

#include <iostream>

using namespace gemmi;
using namespace RDKit;

int main(int argc, char const *argv[]) {
  auto startTotTime = std::chrono::high_resolution_clock::now();
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <cif file>" << std::endl;
    return 1;
  }
  std::string file_name = argv[1];
  bool flag_b = false;

  for (int i = 2; i < argc; ++i) {
    if (std::string(argv[i]) == "-b") {
      flag_b = true;
    }
  }

  Structure st = read_structure_gz(file_name);
  std::cout << "Done: Loading structure" << std::endl;

  RDKit::Conformer *conf = new RDKit::Conformer();
  RDKit::RWMol mol = gemmiStructureToRDKit(st, *conf, false);
  mol.addConformer(conf, true);

  double cutoff = 4.5; // 01;
  std::vector<RDGeom::Point3D> atom_coords = conf->getPositions();
  FastNS grid(atom_coords, cutoff);
  auto results = grid.selfSearch();
  std::cout << "Done: Computing Neighbors" << std::endl;

  if (!flag_b) {
    perceiveBonds(mol, results, 0.45);

    // NOTE: should we clean up bonds before or after adding the struct conns?
    mol.updatePropertyCache(false);
    // Clean up incorrect bonds
    CleanUpMolecule(mol, *conf);

    for (gemmi::Connection &conn : st.connections) {
      // FIX: Need to iterate over all models
      gemmi::Atom *a1 = st.first_model().find_cra(conn.partner1).atom;
      gemmi::Atom *a2 = st.first_model().find_cra(conn.partner2).atom;

      if (mol.getBondBetweenAtoms(a1->serial - 1, a2->serial - 1) == nullptr) {
        mol.addBond((unsigned int)a1->serial - 1, (unsigned int)a2->serial - 1,
                    RDKit::Bond::BondType::SINGLE);
        continue;
      }
    }

    for (auto atomIt = mol.beginAtoms(); atomIt != mol.endAtoms(); ++atomIt) {
      auto atom = *atomIt;
      atom->setHybridization(HybridizationType::SP);
    }
    PerceiveBondOrders(mol);
  }
  // for (gemmi::Connection &conn : st.connections) {
  //   // FIX: Need to iterate over all models
  //   gemmi::Atom *a1 = st.first_model().find_cra(conn.partner1).atom;
  //   gemmi::Atom *a2 = st.first_model().find_cra(conn.partner2).atom;
  //
  //   if (mol.getBondBetweenAtoms(a1->serial - 1, a2->serial - 1) == nullptr) {
  //     mol.addBond((unsigned int)a1->serial - 1, (unsigned int)a2->serial - 1,
  //                 RDKit::Bond::BondType::SINGLE);
  //     continue;
  //   }
  // }
  //

  if (flag_b) {
    findBondsDeconstructedRDKit(mol, results);
    std::cout << "Done: Finding bond orders" << std::endl;
  }

  std::string bondOrders = "";
  for (auto bondIt = mol.beginBonds(); bondIt != mol.endBonds(); ++bondIt) {
    RDKit::Bond *bond = *bondIt;
    auto res1 = mol.getAtomWithIdx(bond->getBeginAtomIdx())->getMonomerInfo();
    auto res2 = mol.getAtomWithIdx(bond->getEndAtomIdx())->getMonomerInfo();

    AtomPDBResidueInfo *residue1 = dynamic_cast<AtomPDBResidueInfo *>(res1);
    AtomPDBResidueInfo *residue2 = dynamic_cast<AtomPDBResidueInfo *>(res2);

    bondOrders += std::to_string(bond->getBeginAtomIdx()) + " " +
                  // res1->getName() + " " + residue1->getResidueName() + " " +
                  // residue1->getAltLoc() + " " +
                  std::to_string(bond->getEndAtomIdx()) + " " +
                  // res2->getName() + " " + residue2->getResidueName() + " " +
                  // residue2->getAltLoc() + " " +
                  std::to_string(bond->getBondType()) + "\n";
  }
  std::cout << "FINAL RESULT: " << std::endl;
  std::cout << bondOrders << std::endl;

  return 0;
}
