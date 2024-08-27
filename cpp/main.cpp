#include "lahuta.hpp"
#include <gemmi/mmread_gz.hpp> // for read_structure_gz

#include <chrono>
#include <iostream>
#include <string>

#include "test.hpp"

#define T() std::chrono::high_resolution_clock::now()
#define TO_MS(d) std::chrono::duration_cast<std::chrono::milliseconds>(d)

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
  lahuta::Luni luni(file_name);

  auto neighbors = luni.find_neighbors();

  RDKit::RWMol *mol = &luni.get_molecule();


  auto atom_iter_start = T();
  for (auto &atom: mol->atoms()) {
    auto *res1 = dynamic_cast<RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
  };
  auto iterTime = TO_MS(T() - atom_iter_start).count();
  std::cout << "Atom ITER TIME: " << iterTime << " ms" << std::endl;



  auto& atoms = getAtoms();

  // Create a generator with a function that prints the atom index
  AtomGenerator<zAtom> gen(atoms, [](zAtom& atom) {
    std::cout << "Atom index: " << atom.getIdx() << std::endl;
  });

  // Iterate over the generator
  for (auto& atom : gen) {
    // Do something with each atom after the function is applied, if needed
    // std::cout << "Atom index: " << atom.getIdx() << std::endl;
  }


  // std::string smarts_test = "[a;r5,!R1&r4,!R1&r3]1:[a;r5,!R1&r4,!R1&r3]:[a;r5,!R1&r4,!R1&r3]:[a;r5,!R1&r4,!R1&r3]:[a;r5,!R1&r4,!R1&r3]:1";
  // auto match = luni.match_smarts_string(smarts_test);
  // // RDKit::RWMol *smarts_mol = RDKit::SmartsToMol(smarts_test);
  // // auto match = RDKit::SubstructMatch(*mol, *smarts_mol);
  // std::cout << "Match: " << match.size() << std::endl;
  // // log match indices
  // for (auto &m : match) {
  //   auto atom = mol->getAtomWithIdx(m[0].second);
  //   std::cout << "Match: " << atom->getIdx() << " " << atom->getSymbol() << std::endl;
  // }

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
