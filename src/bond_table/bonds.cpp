#include "bond_table/bonds.hpp"

int get_intra_bond_order(const std::string &comp_id,
                         const std::string *atom_id1,
                         const std::string *atom_id2) {

  if (*atom_id1 > *atom_id2)
    std::swap(atom_id1, atom_id2);

  if (get_bond_order(comp_id, *atom_id1, *atom_id2))
    return 2;

  if (*atom_id2 == "O" && *atom_id1 == "C") { // O is less common
    if (getResName(comp_id))                  // comp_id is amino acid
      return 2;
  }

  // FIX: test if gperf optimization is needed
  if (combined_base_names.find(comp_id) != combined_base_names.end() &&
      *atom_id1 == "OP1" && *atom_id2 == "P")
    return 2;

  // Default bond order
  return 1;
}
