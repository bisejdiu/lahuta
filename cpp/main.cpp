#include "contacts/interactions.hpp"
#include "file_system.hpp"
#include "lahuta.hpp"

using namespace lahuta;

int main(int argc, char const *argv[]) {
  /*if (argc < 2) {*/
  /*  std::cerr << "Usage: " << argv[0] << " <cif file>" << std::endl;*/
  /*  return 1;*/
  /*}*/

  auto start = std::chrono::high_resolution_clock::now();
  std::string file_name = argv[1];
  /*std::string file_name = "/Users/bsejdiu/projects/lahuta/cpp/data/1kx2_small.cif";*/
  Luni luni(file_name);
  /*luni.assign_arpeggio_atom_types();*/
  /*luni.assign_molstar_atom_types();*/

  if (!luni.success) {
    std::cerr << "Failed to process file: " << file_name << std::endl;
    return 1;
  }
  auto mol = &luni.get_molecule();
  std::cout << "Molecule: " << mol->getNumAtoms() << std::endl;

  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  std::cout << "Time: " << duration.count() << " ms" << std::endl;
  /*return 0;*/

  /*luni.assign_molstar_atom_types();*/

  InteractionOptions opts{5.0};
  Interactions interactions(luni, opts);
  
  std::cout << "HBonds" << std::endl;
  auto _1 = interactions.hbond();
  /*_1.sort_interactions();*/
  /*_1.print_interactions();*/
  std::cout << "size: HBonds: " << _1.size() << std::endl;

  std::cout << "Weak HBonds" << std::endl;
  auto _2 = interactions.weak_hbond();
  /*_2.sort_interactions();*/
  /*_2.print_interactions();*/
  std::cout << "size: Weak HBonds: " << _2.size() << std::endl;

  std::cout << "Hydrophobic" << std::endl;
  auto _3 = interactions.hydrophobic();
  /*_3.sort_interactions();*/
  /*_3.print_interactions();*/
  std::cout << "size: Hydrophobic: " << _3.size() << std::endl;

  std::cout << "Halogen" << std::endl;
  auto _4 = interactions.halogen();
  /*_4.sort_interactions();*/
  /*_4.print_interactions();*/
  std::cout << "Halogen: " << _4.size() << std::endl;

  std::cout << "Ionic" << std::endl;
  auto _5 = interactions.ionic();
  /*_5.sort_interactions();*/
  /*_5.print_interactions();*/
  std::cout << "size: Ionic: " << _5.size() << std::endl;

  std::cout << "Metalic" << std::endl;
  auto _6 = interactions.metalic();
  /*_6.sort_interactions();*/
  /*_6.print_interactions();*/
  std::cout << "size: Metalic: " << _6.size() << std::endl;

  std::cout << "CationPi" << std::endl;
  auto _7 = interactions.cationpi();
  /*_7.sort_interactions();*/
  /*_7.print_interactions();*/
  std::cout << "size: CationPi: " << _7.size() << std::endl;

  std::cout << "PiStacking" << std::endl;
  auto _8 = interactions.pistacking();
  /*_8.sort_interactions();*/
  /*_8.print_interactions();*/
  std::cout << "size: PiStacking: " << _8.size() << std::endl;

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
