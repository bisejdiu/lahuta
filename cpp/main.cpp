#include "lahuta.hpp"
#include "logging.hpp"
#include "selections/tokenizer.hpp"
#include <GraphMol/BondIterators.h>
#include <iostream>
#include <ostream>
#include "contacts/engine.hpp"
#include "contacts/molstar/provider.hpp"
#include "entities/formatter.hpp"

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

  // luni.set_atom_typing_method(ContactComputerType::Arpeggio);
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

  auto &topology = luni.get_topology();
  InteractionEngine<MolStarContactProvider> engine;

  // Compute and log all interaction types
  std::cout << "Computing HBonds..." << std::endl;
  auto start_hb = std::chrono::high_resolution_clock::now();
  auto hbonds = engine.compute(topology, InteractionType::HydrogenBond);
  auto end_hb = std::chrono::high_resolution_clock::now();
  auto hb_duration = std::chrono::duration_cast<std::chrono::microseconds>(end_hb - start_hb);
  std::cout << "HBonds computation time: " << hb_duration.count() << " us" << std::endl;
  ContactTableFormatter::print_contact_table(hbonds, topology, "HydrogenBond");

  std::cout << "Computing Weak HBonds..." << std::endl;
  auto start_weak_hb = std::chrono::high_resolution_clock::now();
  auto weak_hbonds = engine.compute(topology, InteractionType::WeakHydrogenBond);
  auto end_weak_hb = std::chrono::high_resolution_clock::now();
  auto weak_hb_duration = std::chrono::duration_cast<std::chrono::microseconds>(end_weak_hb - start_weak_hb);
  std::cout << "Weak HBonds computation time: " << weak_hb_duration.count() << " us" << std::endl;
  ContactTableFormatter::print_contact_table(weak_hbonds, topology, "WeakHydrogenBond");

  std::cout << "Computing Hydrophobic..." << std::endl;
  auto hydrophobic_start = std::chrono::high_resolution_clock::now();
  auto hydrophobic = engine.compute(topology, InteractionType::Hydrophobic);
  auto hydrophobic_end = std::chrono::high_resolution_clock::now();
  auto hydrophobic_duration = std::chrono::duration_cast<std::chrono::microseconds>(hydrophobic_end - hydrophobic_start);
  std::cout << "Hydrophobic computation time: " << hydrophobic_duration.count() << " us" << std::endl;
  ContactTableFormatter::print_contact_table(hydrophobic, topology, "Hydrophobic");

  std::cout << "Computing Halogen..." << std::endl;
  auto start_halogen = std::chrono::high_resolution_clock::now();
  auto halogen = engine.compute(topology, InteractionType::Halogen);
  auto end_halogen = std::chrono::high_resolution_clock::now();
  auto halogen_duration = std::chrono::duration_cast<std::chrono::microseconds>(end_halogen - start_halogen);
  std::cout << "Halogen computation time: " << halogen_duration.count() << " us" << std::endl;
  ContactTableFormatter::print_contact_table(halogen, topology, "Halogen");

  std::cout << "Computing Ionic..." << std::endl;
  auto start_ionic = std::chrono::high_resolution_clock::now();
  auto ionic = engine.compute(topology, InteractionType::Ionic);
  auto end_ionic = std::chrono::high_resolution_clock::now();
  auto ionic_duration = std::chrono::duration_cast<std::chrono::microseconds>(end_ionic - start_ionic);
  std::cout << "Ionic computation time: " << ionic_duration.count() << " us" << std::endl;
  ContactTableFormatter::print_contact_table(ionic, topology, "Ionic");

  std::cout << "Computing MetalCoordination..." << std::endl;
  auto start_metalic = std::chrono::high_resolution_clock::now();
  auto metalic = engine.compute(topology, InteractionType::MetalCoordination);
  auto end_metalic = std::chrono::high_resolution_clock::now();
  auto metalic_duration = std::chrono::duration_cast<std::chrono::microseconds>(end_metalic - start_metalic);
  std::cout << "MetalCoordination computation time: " << metalic_duration.count() << " us" << std::endl;
  ContactTableFormatter::print_contact_table(metalic, topology, "MetalCoordination");

  std::cout << "Computing CationPi..." << std::endl;
  auto start_cationpi = std::chrono::high_resolution_clock::now();
  auto cationpi = engine.compute(topology, InteractionType::CationPi);
  auto end_cationpi = std::chrono::high_resolution_clock::now();
  auto cationpi_duration = std::chrono::duration_cast<std::chrono::microseconds>(end_cationpi - start_cationpi);
  std::cout << "CationPi computation time: " << cationpi_duration.count() << " us" << std::endl;
  ContactTableFormatter::print_contact_table(cationpi, topology, "CationPi");

  std::cout << "Computing PiStacking..." << std::endl;
  auto start_pistacking = std::chrono::high_resolution_clock::now();
  auto pistacking = engine.compute(topology, InteractionType::PiStacking);
  auto end_pistacking = std::chrono::high_resolution_clock::now();
  auto pistacking_duration = std::chrono::duration_cast<std::chrono::microseconds>(end_pistacking - start_pistacking);
  std::cout << "PiStacking computation time: " << pistacking_duration.count() << " us" << std::endl;
  ContactTableFormatter::print_contact_table(pistacking, topology, "PiStacking");

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
