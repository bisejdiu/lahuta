#include "lahuta.hpp"
#include <gemmi/mmread_gz.hpp> // for read_structure_gz

#include <chrono>
#include <iostream>
#include <string>

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

  // auto load_start = T();
  // Structure st = read_structure_gz(file_name);
  // auto loadTime = TO_MS(T() - load_start).count();
  // std::cout << "Load: " << loadTime << "ms" << std::endl;

  // Current API
  // auto source = Lahuta::GemmiSource();
  // source.process(file_name);

  // Lahuta::Luni luni(source);
  Lahuta::Luni luni(file_name);

  // FIX: replace getNeighborPairsSize with just Size()
  std::cout << "Neighbors: " << luni.get_neighbors().size() << std::endl;

  auto neighbors = luni.get_neighbors();
  RDKit::RWMol *mol = &luni.get_molecule();

  auto log_bond_info = [&](const RDKit::Bond *bond) {
    auto first_atom = mol->getAtomWithIdx(bond->getBeginAtomIdx());
    auto second_atom = mol->getAtomWithIdx(bond->getEndAtomIdx());

    // residue info
    auto *res1 =
        dynamic_cast<RDKit::AtomPDBResidueInfo *>(first_atom->getMonomerInfo());
    auto *res2 = dynamic_cast<RDKit::AtomPDBResidueInfo *>(
        second_atom->getMonomerInfo());

    auto atom1_name = res1->getName();
    auto atom2_name = res2->getName();
    //
    std::string residue1 =
        res1->getResidueName() + "-" + std::to_string(res1->getResidueNumber());
    std::string residue2 =
        res2->getResidueName() + "-" + std::to_string(res2->getResidueNumber());

    auto bond_order = std::to_string(bond->getBondTypeAsDouble());
    if (bond->getIsAromatic()) {
      bond_order = "a";
    }

    std::cout << " " << residue1 << " " << residue2 << " "
              << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx() << " "
              << atom1_name << " " << atom2_name << " "
              << bond_order << std::endl;
  };

  int o1 = 0;
  int o2 = 0;
  int aromatic = 0;
  for (auto bondIt = mol->beginBonds(); bondIt != mol->endBonds(); ++bondIt) {
    RDKit::Bond *bond = *bondIt;

    if (bond->getBondType() == RDKit::Bond::BondType::SINGLE) {
      // log_bond_info(bond);
      o1++;
    } else if (bond->getBondType() == RDKit::Bond::BondType::DOUBLE) {
      // log_bond_info(bond);
      o2++;
    }
    if (bond->getIsAromatic()) {
      // log_bond_info(bond);
      aromatic++;
    }
  }
  std::cout << "1: " << o1 << " 2: " << o2 << " a: " << aromatic
            << " t: " << o1 + o2 + aromatic << std::endl;
  std::cout << "Nr. Bonds: " << mol->getNumBonds() << std::endl;

  auto totTime = TO_MS(T() - startTotTime).count();
  std::cout << "Total Time: " << totTime << "ms" << std::endl;

  return 0;
}
