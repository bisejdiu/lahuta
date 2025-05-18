#include "contacts/interactions.hpp"
#include "lahuta.hpp"
#include "logging.hpp"
#include "selections/tokenizer.hpp"
#include <GraphMol/BondIterators.h>

using namespace lahuta;

int main(int argc, char const *argv[]) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <cif file>" << std::endl;
    return 1;
  }

  Logger::get_instance().set_log_level(Logger::LogLevel::Trace);

  auto start_t = std::chrono::high_resolution_clock::now();

   std::string file_name = argv[1];

  auto start = std::chrono::high_resolution_clock::now();
  Luni luni(file_name);
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
  

  InteractionOptions opts{5.0};
  auto &topology = luni.get_topology();
  Interactions interactions(topology, opts);

  std::cout << "HBonds" << std::endl;
  auto start_hb = std::chrono::high_resolution_clock::now();
  auto _1 = interactions.hbond();
  auto end_hb = std::chrono::high_resolution_clock::now();
  auto hb_duration = std::chrono::duration_cast<std::chrono::microseconds>(end_hb - start_hb);
  std::cout << "Time: " << hb_duration.count() << " us" << std::endl;
  /*_1.sort_interactions();*/
  /*_1.print_interactions();*/
  std::cout << "size: HBonds: " << _1.size() << std::endl;

  std::cout << "Weak HBonds" << std::endl;
  auto start_weak_hb = std::chrono::high_resolution_clock::now();
  auto _2 = interactions.weak_hbond();
  auto end_weak_hb = std::chrono::high_resolution_clock::now();
  auto weak_hb_duration = std::chrono::duration_cast<std::chrono::microseconds>(end_weak_hb - start_weak_hb);
  std::cout << "Time: " << weak_hb_duration.count() << " us" << std::endl;
  /*_2.sort_interactions();*/
  /*_2.print_interactions();*/
  std::cout << "size: Weak HBonds: " << _2.size() << std::endl;

  std::cout << "Hydrophobic" << std::endl;
  auto hb_start = std::chrono::high_resolution_clock::now();
  auto _3 = interactions.hydrophobic();
  auto hb_end = std::chrono::high_resolution_clock::now();
  auto hb_duration2 = std::chrono::duration_cast<std::chrono::microseconds>(hb_end - hb_start);
  std::cout << "Time: " << hb_duration2.count() << " us" << std::endl;
  /*_3.sort_interactions();*/
  /*_3.print_interactions();*/
  std::cout << "size: Hydrophobic: " << _3.size() << std::endl;

  std::cout << "Halogen" << std::endl;
  auto start_halogen = std::chrono::high_resolution_clock::now();
  auto _4 = interactions.halogen();
  auto end_halogen = std::chrono::high_resolution_clock::now();
  auto halogen_duration = std::chrono::duration_cast<std::chrono::microseconds>(end_halogen - start_halogen);
  std::cout << "Time: " << halogen_duration.count() << " us" << std::endl;
  /*_4.sort_interactions();*/
  /*_4.print_interactions();*/
  // std::cout << "size: Halogen: " << _4.size() << std::endl;

  std::cout << "Ionic" << std::endl;
  auto start_ionic = std::chrono::high_resolution_clock::now();
  auto _5 = interactions.ionic();
  auto end_ionic = std::chrono::high_resolution_clock::now();
  auto ionic_duration = std::chrono::duration_cast<std::chrono::microseconds>(end_ionic - start_ionic);
  std::cout << "Time: " << ionic_duration.count() << " us" << std::endl;
  /*_5.sort_interactions();*/
  /*_5.print_interactions();*/
  std::cout << "size: Ionic: " << _5.size() << std::endl;

  std::cout << "Metalic" << std::endl;
  auto start_metalic = std::chrono::high_resolution_clock::now();
  auto _6 = interactions.metal();
  auto end_metalic = std::chrono::high_resolution_clock::now();
  auto metalic_duration = std::chrono::duration_cast<std::chrono::microseconds>(end_metalic - start_metalic);
  std::cout << "Time: " << metalic_duration.count() << " us" << std::endl;
  /*_6.sort_interactions();*/
  /*_6.print_interactions();*/
  std::cout << "size: Metalic: " << _6.size() << std::endl;

  std::cout << "CationPi" << std::endl;
  auto start_cationpi = std::chrono::high_resolution_clock::now();
  auto _7 = interactions.cationpi();
  auto end_cationpi = std::chrono::high_resolution_clock::now();
  auto cationpi_duration = std::chrono::duration_cast<std::chrono::microseconds>(end_cationpi - start_cationpi);
  std::cout << "Time: " << cationpi_duration.count() << " us" << std::endl;
  /*_7.sort_interactions();*/
  /*_7.print_interactions();*/
  std::cout << "size: CationPi: " << _7.size() << std::endl;

  std::cout << "PiStacking" << std::endl;
  auto start_pistacking = std::chrono::high_resolution_clock::now();
  auto _8 = interactions.pistacking();
  auto end_pistacking = std::chrono::high_resolution_clock::now();
  auto pistacking_duration = std::chrono::duration_cast<std::chrono::microseconds>(end_pistacking - start_pistacking);
  std::cout << "Time: " << pistacking_duration.count() << " us" << std::endl;
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
