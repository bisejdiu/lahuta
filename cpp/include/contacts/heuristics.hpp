/*
    This code adapts functionality from the Mol* (mol-star) molecular visualization package,
    originally in TypeScript, for use in C++ with RDKit in the Lahuta project.

    Mol* is licensed under the MIT License (https://github.com/molstar/molstar),
    while this adaptation is released under the GNU General Public License (GPL).
*/
#ifndef LAHUTA_HEURISTICS_HPP
#define LAHUTA_HEURISTICS_HPP

#include "utils.hpp"
#include <GraphMol/Atom.h>
#include <GraphMol/RWMol.h>

namespace lahuta {

/// Carbon in a carboxylate group
inline bool is_carboxylate(const RDKit::RWMol &mol, const RDKit::Atom &atom) {

  if (atom.getAtomicNum() != 6) return false;

  // carbon bonded to exactly two oxygens and one carbon
  if ((get_bond_count(mol, atom, 6) != 1) || (get_bond_count(mol, atom, 8) != 2)) return false;

  // no. of terminal oxygens (those with exactly one non-hydrogen bond)
  unsigned int terminal_oxygen_count = 0;
  for (const auto &bond : mol.atomBonds(&atom)) {
    const RDKit::Atom *neighbor_atom = bond->getOtherAtom(&atom);

    if (neighbor_atom->getAtomicNum() == 8) {
      auto bonds_to_non_H = get_bond_count(mol, *neighbor_atom) - get_bond_count(mol, *neighbor_atom, 1);

      if (bonds_to_non_H == 1) {
        terminal_oxygen_count++;
      }
    }
  }

  // two terminal oxygens
  return terminal_oxygen_count == 2;
}

/// Sulfur in a sulfonic acid or sulfonate group
inline bool is_sulfonic_acid(const RDKit::RWMol &mol, const RDKit::Atom &atom) {
  if (atom.getAtomicNum() != 16) return false;

  // bonded to exactly three oxygens
  return get_bond_count(mol, atom, 8) == 3;
}

/// Sulfur in a sulfate group
inline bool is_sulfate(const RDKit::RWMol &mol, const RDKit::Atom &atom) {
  if (atom.getAtomicNum() != 16) return false;

  // bonded to exactly four oxygens?
  return get_bond_count(mol, atom, 8) == 4;
}

/// Phosphorus in a phosphate group
inline bool is_phosphate(const RDKit::RWMol &mol, const RDKit::Atom &atom) {
  if (atom.getAtomicNum() != 15) return false;

  // no. of bonds to the phosphorus atom
  unsigned int total_bonds = get_bond_count(mol, atom);

  // no. of bonds to oxygen
  unsigned int oxygen_bonds = get_bond_count(mol, atom, 8);

  // all bonds are with oxygen
  return oxygen_bonds == total_bonds;
}

/// Carbon in a guanidine group
inline bool is_guanidine(const RDKit::RWMol &mol, const RDKit::Atom &atom) {
  if (atom.getAtomicNum() != 6) return false;

  // three bonds and all to nitrogen
  if (get_bond_count(mol, atom) != 3 || get_bond_count(mol, atom, 7) != 3) return false;

  // no. of terminal nitrogens (exactly one non-hydrogen bond)
  unsigned int terminal_nitrogen_count = 0;
  for (const auto &bond : mol.atomBonds(&atom)) {
    const RDKit::Atom *neighbor_atom = bond->getOtherAtom(&atom);

    if (neighbor_atom->getAtomicNum() == 7) {
      int bonds_to_non_H = get_bond_count(mol, *neighbor_atom) - get_bond_count(mol, *neighbor_atom, 1);

      if (bonds_to_non_H == 1) {
        terminal_nitrogen_count++;
      }
    }
  }

  // two terminal nitrogens
  return terminal_nitrogen_count == 2;
}

/// Carbon in an acetamidine group
inline bool is_acetamidine(const RDKit::RWMol &mol, const RDKit::Atom &atom) {
  if (atom.getAtomicNum() != 6) return false;

  // not the most efficient, but reads well
  if (get_bond_count(mol, atom) != 3 ||    // 3 bonds
      get_bond_count(mol, atom, 7) != 2 || // 2 to nitrogen
      get_bond_count(mol, atom, 6) != 1) { // 1 to carbon
    return false;
  }

  // no. of terminal nitrogens (exactly one non-hydrogen bond)
  unsigned int terminal_nitrogen_count = 0;
  for (const auto &bond : mol.atomBonds(&atom)) {
    const RDKit::Atom *neighbor_atom = bond->getOtherAtom(&atom);

    if (neighbor_atom->getAtomicNum() == 7) {
      auto bonds_to_non_H = get_bond_count(mol, *neighbor_atom) - get_bond_count(mol, *neighbor_atom, 1);

      if (bonds_to_non_H == 1) {
        terminal_nitrogen_count++;
      }
    }
  }

  // two terminal nitrogens
  return terminal_nitrogen_count == 2;
}

} // namespace lahuta

#endif // LAHUTA_HEURISTICS_HPP
