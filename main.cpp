#include <gemmi/mmread_gz.hpp> // for read_structure_gz
#include "lahuta.hpp"

#include <chrono>
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

  auto prog_start = T();
  auto load_start = T();
  Structure st = read_structure_gz(file_name);
  auto load_end = T();
  std::cout << "Done: Loading structure" << std::endl;
  std::cout << "perf: " << TO_MS(load_end - load_start).count() << "ms"
            << std::endl;

  Lahuta lahuta(st);
  RDKit::RWMol mol = lahuta.getMol();

  // for (gemmi::Connection &conn : st.connections) {
  //   // FIX: Need to iterate over all models
  //   gemmi::Atom *a1 = st.first_model().find_cra(conn.partner1).atom;
  //   gemmi::Atom *a2 = st.first_model().find_cra(conn.partner2).atom;
  //
  //   if (mol.getBondBetweenAtoms(a1->serial - 1, a2->serial - 1) == nullptr) {
  //     // compute distance between the two atoms:
  //     double dist =
  //         (atom_coords[a1->serial - 1] - atom_coords[a2->serial - 1]).length();
  //     std::cout << "conn pair: " << a1->serial << " " << a2->serial << " "
  //               << dist << std::endl;
  //     // mol.addBond((unsigned int)a1->serial - 1, (unsigned int)a2->serial-1,
  //     //             RDKit::Bond::BondType::SINGLE);
  //     // continue;
  //   }
  // }

  // std::string bO = "";
  int FirstOrderCount = 0;
  int SecondOrderCount = 0;
  for (auto bondIt = mol.beginBonds(); bondIt != mol.endBonds(); ++bondIt) {
    RDKit::Bond *bond = *bondIt;

    if (bond->getBondType() == RDKit::Bond::BondType::SINGLE) {
      FirstOrderCount++;
    } else if (bond->getBondType() == RDKit::Bond::BondType::DOUBLE) {
      SecondOrderCount++;
    }

    //   auto res1 =
    //   mol.getAtomWithIdx(bond->getBeginAtomIdx())->getMonomerInfo(); auto
    //   res2 = mol.getAtomWithIdx(bond->getEndAtomIdx())->getMonomerInfo();
    //
    //   AtomPDBResidueInfo *residue1 = dynamic_cast<AtomPDBResidueInfo
    //   *>(res1); AtomPDBResidueInfo *residue2 =
    //   dynamic_cast<AtomPDBResidueInfo *>(res2);
    //
    //   bO += std::to_string(bond->getBeginAtomIdx()) + " " +
    //         // res1->getName() + " " + residue1->getResidueName() + " " +
    //         // residue1->getAltLoc() + " " +
    //         std::to_string(bond->getEndAtomIdx()) + " " +
    //         // res2->getName() + " " + residue2->getResidueName() + " " +
    //         // residue2->getAltLoc() + " " +
    //         std::to_string(bond->getBondType()) + "\n";
  }
  // std::cout << "FINAL RESULT: " << std::endl;
  // std::cout << bO << std::endl;
  std::cout << "First Order Count: " << FirstOrderCount << std::endl;
  std::cout << "Second Order Count: " << SecondOrderCount << std::endl;
  std::cout << "Total Bonds: " << mol.getNumBonds() << std::endl;

  auto prog_end = T();
  std::cout << "Total Time: " << TO_MS(prog_end - prog_start).count() << "ms"
            << std::endl;

  return 0;
}
