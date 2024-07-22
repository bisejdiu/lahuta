#include "bonds.hpp"
#include "GraphMol/Atom.h"
#include "GraphMol/MonomerInfo.h"

using HybridizationType = RDKit::Atom::HybridizationType;

int get_intra_bond_order(const std::string &comp_id, const std::string *atom_id1,
                         const std::string *atom_id2) {

  // const std::string *atom_id1 = &atom_id1;
  // const std::string *atom_id2 = &atom_id2;

  if (*atom_id1 > *atom_id2)
    std::swap(atom_id1, atom_id2);

  // Use tuple directly to avoid string concatenation
  auto key = std::make_tuple(comp_id, *atom_id1, *atom_id2);
  auto it = intra_bond_order_table.find(key);
  if (it != intra_bond_order_table.end())
    return it->second;

  if (combined_amino_acid_names.find(comp_id) != combined_amino_acid_names.end() && *atom_id1 == "C" && *atom_id2 == "O")
    return 2;
  if (combined_base_names.find(comp_id) != combined_base_names.end() && *atom_id1 == "OP1" && *atom_id2 == "P")
    return 2;

  // Default bond order
  return 1;
}

