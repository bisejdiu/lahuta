#include <gemmi/mmread_gz.hpp> // for read_structure_gz

#include "nsgrid.hpp" // for FastNS
#include "bonds.hpp"  // for perceiveBonds
#include "conv.hpp" // for gemmiStructureToRDKit

#include <chrono>

#include "bond_utils.hpp" // for CleanUpMolecule
#include "bond_order.hpp"

#include <iostream>

using namespace gemmi;
using namespace RDKit;

int main(int argc, char const *argv[]) {
  auto startTotTime = std::chrono::high_resolution_clock::now();
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <cif file>" << std::endl;
    return 1;
  }
  std::string file_name = argv[1];
  auto start = std::chrono::high_resolution_clock::now();
  Structure st = read_structure_gz(file_name);
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::cout << "Time to load structure: " << elapsed.count() * 1000 << " ms"
            << std::endl;

  start = std::chrono::high_resolution_clock::now();
  RDKit::Conformer *conf = new RDKit::Conformer();
  RDKit::RWMol mol = gemmiStructureToRDKit(st, *conf, false);
  mol.addConformer(conf, true);
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
  std::cout << "Time to convert to RDKit: " << elapsed.count() * 1000 << " ms"
            << std::endl;


  double cutoff = 4.5; // 01;
  std::vector<RDGeom::Point3D> atom_coords = conf->getPositions();
  FastNS grid(atom_coords, cutoff);
  auto results = grid.selfSearch();
  std::cout << "Number of neighbors: " << results.getNeighborPairsSize()
            << std::endl;

  // perceives bonds
  start = std::chrono::high_resolution_clock::now();
  perceiveBonds(mol, results, 0.45);
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
  std::cout << "Time to perceive bonds: " << elapsed.count() * 1000 << " ms"
            << std::endl;
  std::cout << "Number of bonds found using NSGrid: " << mol.getNumBonds()
            << std::endl;


  // NOTE: should we clean up bonds before or after adding the struct conns?
  mol.updatePropertyCache(false);
  // Clean up incorrect bonds
  start = std::chrono::high_resolution_clock::now();
  CleanUpMolecule(mol, *conf);
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
  std::cout << "Time to clean up molecule: " << elapsed.count() * 1000 << "ms "
            << std::endl;

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

  std::cout << "Number of bonds w/ _struct_conns " << mol.getNumBonds()
            << std::endl;

  // get ROMol from RWMol using cast
  // auto _mol = static_cast<RDKit::ROMol*>(new RDKit::RWMol(mol));

  for (auto atomIt = mol.beginAtoms(); atomIt != mol.endAtoms(); ++atomIt) {
    auto atom = *atomIt;
    atom->setHybridization(HybridizationType::SP);
  }

  start = std::chrono::high_resolution_clock::now();
  PerceiveBondOrders(mol);
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
  std::cout << "Time to perceive bond orders: " << elapsed.count() * 1000
            << " ms" << std::endl;

  // for (auto atom : mol.atoms()) {
  //   auto resi = atom->getMonomerInfo();
  //   AtomPDBResidueInfo *residue = dynamic_cast<AtomPDBResidueInfo *>(resi);
  //   std::cout << "Atom Info: " << atom->getIdx() << " " << atom->getSymbol()
  //             << " " << resi->getMonomerType() << " " << resi->getName() << "
  //             "
  //             << residue->getResidueName() << " " <<
  //             residue->getIsHeteroAtom()
  //             << " " << std::endl;
  // }

  // std::string bondOrders = "";
  // for (auto bondIt = mol.beginBonds(); bondIt != mol.endBonds(); ++bondIt) {
  //   RDKit::Bond *bond = *bondIt;
  //   bondOrders += std::to_string(bond->getBeginAtomIdx()) + " " +
  //                 std::to_string(bond->getEndAtomIdx()) + " " +
  //                 std::to_string(bond->getBondType()) + "\n";
  // }
  // std::cout << "FINAL RESULT: " << std::endl;
  // std::cout << bondOrders << std::endl;

  auto endTotTime = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsedTotTime = endTotTime - startTotTime;
  std::cout << "Total time: " << elapsedTotTime.count() * 1000 << " ms"
            << std::endl;
  return 0;
}
