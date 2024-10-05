#include <chrono>
#include <iostream>
#include <random>
#include <string>

#include "lahuta.hpp"
#include "neighbors.hpp"
#include "visitor.hpp"

using namespace gemmi;
using namespace RDKit;
using namespace lahuta;

// Define test selection strings
std::vector<std::string> selections = {
    // Test single term with a single value
    "resid3",

    // Test single term with a range
    "resid3 -10",

    // Test 'not' operator with a single term
    "notresid3",

    // Test term with a list of numeric values
    "resid 3 5 7",

    // Test term with a list of string values
    "resname ALA GLY",

    // Test 'and' operator with range and term
    "resid 1 - 50 and resname ASP",

    // Test 'or' operator with range and term
    "resid 30 - 40 or resname LEU",

    // Test 'not' operator with a term
    "not resname ASP",

    // Test 'and' and 'not' operators combined
    "resid10-20andnotresnameLEU",
    /*"resid 10 - 20 and not resname LEU",*/

    // Test 'or' operator with two ranges
    "resid 5 - 15 or resid 20 - 25",

    // Test grouping with parentheses
    "( resid 5 - 10 or resid 20 - 25 ) and resname GLY",

    // Test 'not' operator with grouping
    "not ( resname ALA or resname LEU )",

    // Test 'and' operator with lists
    "resid 1 2 3 and resname GLY",

    // Test 'and' and 'not' operators with lists
    "resid 1 - 100 and not ( resname ASP LEU )",

    // Complex expression with grouping and operators
    "( resid 10 - 20 and resname GLY ) or ( resid 30 - 40 and resname LEU )",

    "resid -3 - 3 and not resid 0 1 2",
    // Test parentheses overriding operator precedence
    // Lahuta parser `and` has higher precedence than `or`
    // MDAnalysis parser will parse by order of appearance (I think?)
    /*"resid 1 - 5 or resid 6 - 10 and resname ALA",*/
    /*"(resid 1 - 5 or resid 6 - 10 ) and resname ALA",*/
    /*"resid 1 - 5 or (resid 6 - 10 and resname ALA)",*/
};

int main(int argc, char const *argv[]) {
  auto start_total_time = std::chrono::high_resolution_clock::now();
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <cif file>" << std::endl;
    return 1;
  }
  std::string file_name = argv[1];

  Luni luni(file_name);

  auto neighbors = luni.find_neighbors(5.0, 1);
  auto vv = neighbors.type_filter(AtomType::AROMATIC, 0);
  auto mol = &luni.get_molecule();

  Pairs p = neighbors.get_pairs();
  Distances d = neighbors.get_distances();

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(0, p.size() - 1);
  std::uniform_real_distribution<float> disf(0.0, 1.0); // random float between 0 and 1

  std::cout << "p size: " << p.size() << "\n";
  std::vector<AtomAtomPair> _p2; // store only the first 20 pairs
  for (size_t i = 0; i < p.size(); ++i) {
    if (i < 20) {
      _p2.push_back(AtomAtomPair(p[i].first, p[i].second, d[i]));
    }
  }
  Neighbors<AtomAtomPair> _n1(neighbors);
  Neighbors<AtomAtomPair> _n2(luni, _p2);
  auto ns = _n1.intersection(_n1.get_data(), _n2.get_data());
  auto _ns_ = _n1.intersection(_n2);
  std::cout << "ns size: " << ns.size() << "\n";
  std::cout << "_ns_ size: " << _ns_.size() << "\n";

  auto nn = luni.find_ring_neighbors(6.0);
  std::cout << "nn size: " << nn.size() << "\n";

  auto log_bond_info = [&](const RDKit::Bond *bond) {
    auto first_atom = mol->getAtomWithIdx(bond->getBeginAtomIdx());
    auto second_atom = mol->getAtomWithIdx(bond->getEndAtomIdx());
    auto *res1 =
        static_cast<RDKit::AtomPDBResidueInfo *>(first_atom->getMonomerInfo());
    auto *res2 =
        static_cast<RDKit::AtomPDBResidueInfo *>(second_atom->getMonomerInfo());
    auto atom1_name = res1->getName();
    auto atom2_name = res2->getName();
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
              << atom1_name << " " << atom2_name << " " << bond_order
              << std::endl;
  };

  int o1{}, o2{}, aromatic{};
  for (auto bond_it = mol->beginBonds(); bond_it != mol->endBonds();
       ++bond_it) {
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
  std::cout << "1: " << o1 << " 2: " << o2 << " a: " << aromatic
            << " t: " << o1 + o2 + aromatic << std::endl;
  std::cout << "Nr. Bonds: " << mol->getNumBonds() << std::endl;

  auto tot_time = to_ms(t() - start_total_time).count();
  std::cout << "Total Time: " << tot_time << "ms" << std::endl;

  return 0;
}
