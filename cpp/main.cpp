#include <chrono>
#include <iostream>
#include <string>

#include "lahuta.hpp"
#include "neighbors.hpp"

using namespace gemmi;
using namespace RDKit;

int main(int argc, char const *argv[]) {
  auto startTotTime = std::chrono::high_resolution_clock::now();
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <cif file>" << std::endl;
    return 1;
  }
  std::string file_name = argv[1];

  lahuta::Luni luni(file_name);
  auto neighbors = luni.find_neighbors<lahuta::AtomAtomPair>(5.0, 1);
  auto mol = &luni.get_molecule();

  std::cout << "Testing filtering" << std::endl;
  std::vector<int> atom_indices;
  for (int i = 0; i < 134; i++) {
    atom_indices.push_back(i);
  }
  lahuta::Luni new_luni = luni.filter_luni(atom_indices);
  auto new_mol = &new_luni.get_molecule();
  std::cout << "Luni n_atoms" << luni.n_atoms() << std::endl; 
  std::cout << "New Mol: " << new_mol->getNumAtoms() << std::endl;
  // first 20
  std::vector<int> atom_indices2 = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 
                                    10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
  lahuta::Luni new_luni2 = new_luni.filter_luni(atom_indices2);
  auto new_mol2 = &new_luni2.get_molecule();
  std::cout << "Luni n_atoms" << luni.n_atoms() << std::endl;
  std::cout << "New Mol2: " << new_mol2->getNumAtoms() << std::endl;
  
  // log aromatic atoms
  for (auto &atom: new_luni.get_molecule().atoms()) {
    if (atom->getIsAromatic()) {
      std::cout << "Aromatic Atom: " << atom->getIdx() << std::endl;
    }
  }
  lahuta::RingDataVec rings = new_luni.get_rings();
  for (auto &ring: rings.rings) {
    std::cout << "Ring: ";
    for (auto &atom: ring.atom_ids) {
      std::cout << atom << " ";
    }
    std::cout << std::endl;
  }



  auto log_bond_info = [&](const RDKit::Bond *bond) {
    auto first_atom = mol->getAtomWithIdx(bond->getBeginAtomIdx());
    auto second_atom = mol->getAtomWithIdx(bond->getEndAtomIdx());
    auto *res1 = static_cast<RDKit::AtomPDBResidueInfo *>(first_atom->getMonomerInfo());
    auto *res2 = static_cast<RDKit::AtomPDBResidueInfo *>(second_atom->getMonomerInfo());
    auto atom1_name = res1->getName();
    auto atom2_name = res2->getName();
    std::string residue1 = res1->getResidueName() + "-" + std::to_string(res1->getResidueNumber());
    std::string residue2 = res2->getResidueName() + "-" + std::to_string(res2->getResidueNumber());

    auto bond_order = std::to_string(bond->getBondTypeAsDouble());
    if (bond->getIsAromatic()) {
      bond_order = "a";
    }

    std::cout << " " << residue1 << " " << residue2 << " "
              << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx() << " "
              << atom1_name << " " << atom2_name << " "
              << bond_order << std::endl;
  };

  int o1{}, o2{}, aromatic{};
  for (auto bondIt = mol->beginBonds(); bondIt != mol->endBonds(); ++bondIt) {
    RDKit::Bond *bond = *bondIt;

    if (bond->getBondType() == RDKit::Bond::BondType::SINGLE) {
      o1++;
    } else if (bond->getBondType() == RDKit::Bond::BondType::DOUBLE) {
      o2++;
    }
    if (bond->getIsAromatic()) {
      aromatic++;
    }
  }
  std::cout << "1: " << o1 << " 2: " << o2 << " a: " << aromatic
            << " t: " << o1 + o2 + aromatic << std::endl;
  std::cout << "Nr. Bonds: " << mol->getNumBonds() << std::endl;

  auto totTime = to_ms(t() - startTotTime).count();
  std::cout << "Total Time: " << totTime << "ms" << std::endl;

  return 0;
}
