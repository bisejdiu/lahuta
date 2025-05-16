#include "contacts/interactions.hpp"
#include "lahuta.hpp"
#include "logging.hpp"
#include "selections/tokenizer.hpp"

using namespace lahuta;

int main(int argc, char const *argv[]) {
  /*if (argc < 2) {*/
  /*  std::cerr << "Usage: " << argv[0] << " <cif file>" << std::endl;*/
  /*  return 1;*/
  /*}*/

  Logger::get_instance().set_log_level(Logger::LogLevel::Trace);

  auto start_t = std::chrono::high_resolution_clock::now();

   std::string file_name = argv[1];

  auto start = std::chrono::high_resolution_clock::now();
  /*std::string file_name = "/Users/bsejdiu/projects/lahuta/cpp/data/1kx2_small.cif";*/
  Luni luni(file_name);
  /*luni.assign_arpeggio_atom_types();*/
  /*luni.assign_molstar_atom_types();*/
  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  std::cout << "Time: " << duration.count() << " us" << std::endl;

  auto start_topology = std::chrono::high_resolution_clock::now();
  if (!luni.build_topology()) {
    std::cerr << "Failed to process file: " << file_name << std::endl;
    return 1;
  }
  auto end_topology = std::chrono::high_resolution_clock::now();
  auto duration_topology = std::chrono::duration_cast<std::chrono::microseconds>(end_topology - start_topology);
  std::cout << "Topology Time: " << duration_topology.count() << " us" << std::endl;

  auto mol = &luni.get_molecule();
  std::cout << "Molecule Atoms: " << mol->getNumAtoms() << std::endl;
  std::cout << "Molecule Bonds: " << mol->getNumBonds() << std::endl;

  // Check if atom typing computation is enabled
  bool is_enabled = luni.is_topology_computation_enabled(TopologyComputation::AtomTyping);
  std::cout << "Atom typing enabled: " << (is_enabled ? "yes" : "no") << std::endl;
  
  // Get atom types using Molstar (already done in topology build)
  const AtomEntityCollection &molstar_types = luni.get_topology().get_atom_types();
  std::cout << "Molstar atom types: " << molstar_types.size() << std::endl;
  
  if (molstar_types.size() > 0) {
    const auto hydrophobic_atoms_molstar = AtomEntityCollection::filter(&luni, AtomType::HYDROPHOBIC);
    std::cout << "Hydrophobic atoms (Molstar): " << hydrophobic_atoms_molstar.size() << std::endl;
  }
  
  // Now try Arpeggio atom typing
  std::cout << "\nSwitching to Arpeggio atom typing..." << std::endl;
  luni.assign_arpeggio_atom_types();
  
  // Check atom types after Arpeggio
  const AtomEntityCollection &arpeggio_types = luni.get_topology().get_atom_types();
  std::cout << "Arpeggio atom types: " << arpeggio_types.size() << std::endl;
  
  if (arpeggio_types.size() > 0) {
    const auto hydrophobic_atoms_arpeggio = AtomEntityCollection::filter(&luni, AtomType::HYDROPHOBIC);
    std::cout << "Hydrophobic atoms (Arpeggio): " << hydrophobic_atoms_arpeggio.size() << std::endl;
  }
  
  // Switch back to Molstar
  std::cout << "\nSwitching back to Molstar atom typing..." << std::endl;
  luni.assign_molstar_atom_types();
  
  // Check atom types after switching back to Molstar
  const AtomEntityCollection &molstar_types_again = luni.get_topology().get_atom_types();
  std::cout << "Molstar atom types (again): " << molstar_types_again.size() << std::endl;
  
  if (molstar_types_again.size() > 0) {
    const auto hydrophobic_atoms_molstar_again = AtomEntityCollection::filter(&luni, AtomType::HYDROPHOBIC);
    std::cout << "Hydrophobic atoms (Molstar again): " << hydrophobic_atoms_molstar_again.size() << std::endl;
  }


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
  for (auto bond: mol->bonds()) {

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
