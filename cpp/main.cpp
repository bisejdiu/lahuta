#include "contacts/interactions.hpp"
#include "file_system.hpp"
#include "lahuta.hpp"

#include "spdlog/common.h"
#include "spdlog/spdlog.h"

using namespace lahuta;

void set_logger_pattern(spdlog::level::level_enum level) {
  if (level == spdlog::level::debug) {
    spdlog::set_pattern("[%T] [%^%l%$] [thread %t] %v");
  } else {
    spdlog::set_pattern("[%^%l%$] %v");
  }
}

int main(int argc, char const *argv[]) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <cif file>" << std::endl;
    return 1;
  }

  auto start = std::chrono::high_resolution_clock::now();
  std::string file_name = argv[1];
  Luni luni(file_name);
  auto mol = &luni.get_molecule();

  if (!luni.success) {
    std::cerr << "Failed to process file: " << file_name << std::endl;
    return 1;
  }
  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  std::cout << "Time: " << duration.count() << " ms" << std::endl;

  /*luni.assign_molstar_atom_types();*/
  InteractionOptions opts{5.0};
  Interactions interactions(luni, opts);
  auto _1 = interactions.hbond();
  _1.sort_interactions();
  _1.print_interactions();
  std::cout << "Weak HBonds" << std::endl;
  auto _2 = interactions.weak_hbond();
  _2.sort_interactions();
  _2.print_interactions();
  auto _3 = interactions.hydrophobic();
  _3.sort_interactions();
  _3.print_interactions();
  auto _4 = interactions.halogen();
  _4.sort_interactions();
  _4.print_interactions();
  auto _5 = interactions.ionic();
  _5.sort_interactions();
  _5.print_interactions();
  auto _6 = interactions.metalic();
  _6.sort_interactions();
  _6.print_interactions();
  auto _7 = interactions.cationpi();
  _7.sort_interactions();
  _7.print_interactions();
  auto _8 = interactions.pistacking();
  _8.sort_interactions();
  _8.print_interactions();

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

    std::cout << " " << residue1 << " " << residue2 << " " << bond->getBeginAtomIdx() << " "
              << bond->getEndAtomIdx() << " " << atom1_name << " " << atom2_name << " " << bond_order
              << std::endl;
  };

  int o1{}, o2{}, aromatic{};
  for (auto bond_it = mol->beginBonds(); bond_it != mol->endBonds(); ++bond_it) {
    RDKit::Bond *bond = *bond_it;

    if (bond->getBondType() == RDKit::Bond::BondType::SINGLE) {
      o1++;
    } else if (bond->getBondType() == RDKit::Bond::BondType::DOUBLE) {
      o2++;
    }
    if (bond->getIsAromatic()) {
      aromatic++;
    }
  }
  std::cout << "1: " << o1 << " 2: " << o2 << " a: " << aromatic << " t: " << o1 + o2 + aromatic << std::endl;
  std::cout << "Nr. Bonds: " << mol->getNumBonds() << std::endl;

  return 0;
}
