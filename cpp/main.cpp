#include <chrono>
#include <iostream>
#include <string>

#include "lahuta.hpp"
#include "neighbors.hpp"
#include "parser.hpp"
#include "visitor.hpp"

using namespace gemmi;
using namespace RDKit;
using namespace lahuta;

//! Parse the expression and return the filtered indices
std::vector<int> parse_and_filter(const Luni &luni,
                                  const std::string &selection) {
  std::vector<std::string> tokens = Luni::tokenize(selection);
  lahuta::Parser parser(tokens);

  // Parse the expression
  lahuta::NodePtr root = parser.parse_expression();

  lahuta::FilterVisitor visitor(luni);
  root->accept(visitor);

  const std::vector<int> &filtered_indices = visitor.get_result();
  return filtered_indices;
}

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
  auto neighbors = luni.find_neighbors<AtomAtomPair>(5.0, 1);
  auto mol = &luni.get_molecule();

  std::cout << "\n\n PARSER TEST \n\n";
  for (const auto &selection : selections) {
    std::cout << "Selection: " << selection << std::endl;

    // Parse and filter the selection
    std::vector<int> indices = parse_and_filter(luni, selection);

    std::cout << "Filtered indices (" << indices.size() << "): ";
    for (int index : indices) {
      std::cout << index << " ";
    }
    std::cout << std::endl << std::endl;
  }

  std::cout << "\n\n PARSER DONE \n\n";

  std::cout << "Testing filtering" << std::endl;
  std::vector<int> atom_indices;
  for (int i = 0; i < 134; i++) {
    atom_indices.push_back(i);
  }
  Luni new_luni = luni.filter_luni(atom_indices);
  auto new_mol = &new_luni.get_molecule();
  std::cout << "Luni n_atoms" << luni.n_atoms() << std::endl;
  std::cout << "New Mol: " << new_mol->getNumAtoms() << std::endl;
  // first 20
  std::vector<int> atom_indices2 = {0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
                                    10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
  Luni new_luni2 = new_luni.filter_luni(atom_indices2);
  auto new_mol2 = &new_luni2.get_molecule();
  std::cout << "Luni n_atoms" << luni.n_atoms() << std::endl;
  std::cout << "New Mol2: " << new_mol2->getNumAtoms() << std::endl;

  // log aromatic atoms
  for (auto &atom : new_luni.get_molecule().atoms()) {
    if (atom->getIsAromatic()) {
      std::cout << "Aromatic Atom: " << atom->getIdx() << std::endl;
    }
  }
  RingDataVec rings = new_luni.get_rings();
  for (auto &ring : rings.rings) {
    std::cout << "Ring: ";
    for (auto &atom : ring.atom_ids) {
      std::cout << atom << " ";
    }
    std::cout << std::endl;
  }

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
