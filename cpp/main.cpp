#include <chrono>
#include <iostream>
#include <string>

#include "GraphMol/RWMol.h"
#include "atom_types.hpp"
#include "contacts/hydrogen_bonds.hpp"
#include "contacts/interactions.hpp"
#include "lahuta.hpp"
#include "neighbors.hpp"
#include "nn.hpp"
#include "nsgrid.hpp"
#include "rings.hpp"
#include "visitor.hpp"

using namespace gemmi;
using namespace RDKit;
using namespace lahuta;

int main(int argc, char const *argv[]) {
  auto start_total_time = std::chrono::high_resolution_clock::now();
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <cif file>" << std::endl;
    return 1;
  }
  std::string file_name = argv[1];

  Luni luni(file_name);

  const auto &G = luni.get_features();

  auto neighbors = luni.find_neighbors(6.0, 0);
  /*auto vv = neighbors.type_filter(AtomType::AROMATIC, 0);*/
  /*auto mol = &luni.get_molecule();*/

  std::cout << "START New contact interface" << std::endl;

  Contacts _c_(&luni);

  auto mol = &luni.get_molecule();
  const auto &atom_entities = luni.get_atom_entities();
  auto atom_neighbors = luni.find_neighbors2(6.0, 10);
  _c_.add_many(atom_neighbors, atom_entities);
  /*_c_.sort_interactions();*/
  /*_c_.print_interactions();*/
  std::cout << "Hbond: " << _c_.size() << std::endl;
  std::cout << "hbond is sorted: " << _c_.is_sorted << std::endl;

  Contacts container(&luni);

  container.add_many(atom_neighbors, atom_entities);
  /*container.sort_interactions();*/
  /*container.print_interactions();*/

  std::vector<RingData> rings = luni.get_rings().rings;
  const auto &ring_entities = luni.get_ring_entities();
  auto ring_neighbors = luni.find_ring_neighbors2(6.0);
  container.add_many(ring_neighbors, ring_entities, atom_entities);

  std::cout << "Atom and Ring Interactions size: " << container.size() << std::endl;
  /*container.sort_interactions();*/
  /*container.print_interactions();*/
  std::cout << "END Atom and Ring Interactions size: \n" << std::endl;

  Contacts hbond_container(&luni);
  GeometryOptions _opts = GeometryOptions();
  NSResults _nb = luni.find_neighbors2(6.0, 0);
  find_hydrogen_bonds(luni, _opts, _nb, hbond_container);
  std::cout << "HBond Interactions size: " << hbond_container.size() << std::endl;
  /*hbond_container.sort_interactions();*/
  /*hbond_container.print_interactions();*/

  Contacts weak_hbond_container(&luni);
  find_weak_hydrogen_bonds(luni, _opts, _nb, weak_hbond_container);
  std::cout << "Weak HBond Interactions size: " << weak_hbond_container.size() << std::endl;
  /*weak_hbond_container.sort_interactions();*/
  /*weak_hbond_container.print_interactions();*/

  /*weak_hbond_container.add(Interaction(0, 7, 300.0, InteractionType::WeakHydrogenBond));*/
  /*weak_hbond_container.sort_interactions();*/
  /*weak_hbond_container.print_interactions();*/

  /*std::cout << "" << std::endl;*/
  /*IC i_sect_1 = hbond_container.set_intersection(weak_hbond_container);*/
  /*i_sect_1.print_interactions();*/
  /*std::cout << "--> Intersection size: " << i_sect_1.size() << std::endl;*/
  /*hbond_container.make_generic();*/
  /*weak_hbond_container.make_generic();*/
  /*IC i_sect_2 = hbond_container.set_intersection(weak_hbond_container);*/
  /*i_sect_2.print_interactions();*/
  /*std::cout << "--> Intersection size: " << i_sect_2.size() << std::endl;*/

  std::cout << "Number of rings: " << rings.size() << std::endl;

  std::cout << "END New contact interface" << std::endl;

  std::vector<AtomType> atypes = luni.get_atom_types();

  std::cout << "STARTING ION CONTACT COMPUTATION" << std::endl;

  std::vector<Feature> group_features = GroupTypeAnalysis::analyze(*mol);
  Interactions processor(&luni, group_features, InteractionOptions{5.0});
  auto ionic = processor.find_ionic_interactions();
  std::cout << "Ionic Interactions size: " << ionic.size() << std::endl;

  /*(ionic & ionic[0]).print_interactions();*/

  /*auto r = ionic[2];*/
  /*std::cout << "ionic &= r: " << (ionic &= r).size() << std::endl;*/
  /*ionic |= ionic[0];*/
  /**/
  /*(ionic |= r).print_interactions();*/
  /*(ionic |= ionic[0]);*/

  container.add(ionic);
  container.sort_interactions();
  container.print_interactions();

  std::cout << "New HBond: " << std::endl;
  auto h2 = processor.find_hbond_interactions();
  /*h2.sort_interactions();*/
  /*h2.print_interactions();*/
  std::cout << "HBond Interactions size: " << h2.size() << std::endl;
  if (h2 == _c_) {
    std::cout << "EQUAL:\n";
  }

  _c_ += h2;
  _c_ ^= h2;
  _c_ |= h2;
  _c_ &= h2;
  _c_ -= h2;

  /*std::cout << "All contacts:" << std::endl;*/
  /*container.print_interactions();*/

  /*Pairs p = neighbors.get_pairs();*/
  /*Distances d = neighbors.get_distances();*/
  /**/
  /*std::random_device rd;*/
  /*std::mt19937 gen(rd());*/
  /*std::uniform_int_distribution<> dis(0, p.size() - 1);*/
  /*std::uniform_real_distribution<float> disf(0.0, 1.0); // random float between 0 and 1*/
  /**/
  /*std::cout << "p size: " << p.size() << "\n";*/
  /*std::vector<AtomAtomPair> _p2; // store only the first 20 pairs*/
  /*for (size_t i = 0; i < p.size(); ++i) {*/
  /*  if (i < 20) {*/
  /*    _p2.push_back(AtomAtomPair(p[i].first, p[i].second, d[i]));*/
  /*  }*/
  /*}*/
  /*Neighbors<AtomAtomPair> _n1(neighbors);*/
  /*Neighbors<AtomAtomPair> _n2(luni, _p2);*/
  /*auto ns = _n1.intersection(_n1.get_data(), _n2.get_data());*/
  /*auto _ns_ = _n1.intersection(_n2);*/
  /*std::cout << "ns size: " << ns.size() << "\n";*/
  /*std::cout << "_ns_ size: " << _ns_.size() << "\n";*/
  /**/
  /*auto nn = luni.find_ring_neighbors(6.0);*/
  /*std::cout << "nn size: " << nn.size() << "\n";*/

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

  auto tot_time = to_ms(t() - start_total_time).count();
  std::cout << "Total Time: " << tot_time << "ms" << std::endl;

  return 0;
}
