#ifndef LAHUTA_HEURISTICS_HPP
#define LAHUTA_HEURISTICS_HPP

#include <GraphMol/Atom.h>
#include <GraphMol/RWMol.h>

#include "chemistry/elements.hpp"
#include "chemistry/utils.hpp"

// clang-format off
namespace lahuta {

inline bool is_C_in_carboxylate(const RDKit::RWMol &mol, const RDKit::Atom &atom) {

  if (atom.getAtomicNum() != Element::C) return false;

  // carbon bonded to exactly two oxygens and one carbon
  if (get_bond_count(mol, atom, Element::C) != 1 || get_bond_count(mol, atom, Element::O) != 2) return false;

  // no. of terminal oxygens (those with exactly one non-hydrogen bond)
  unsigned int terminal_O_count = 0;
  for (const auto &bond : mol.atomBonds(&atom)) {
    const RDKit::Atom *nbr_atom = bond->getOtherAtom(&atom);

    if (nbr_atom->getAtomicNum() == Element::O) {
      auto bonds_to_non_H = get_bond_count(mol, *nbr_atom) - get_bond_count(mol, *nbr_atom, 1);

      if (bonds_to_non_H == 1) {
        terminal_O_count++;
      }
    }
  }

  return terminal_O_count == 2;
}

inline bool is_S_in_sulfonic_acid(const RDKit::RWMol &mol, const RDKit::Atom &atom) {
  if (atom.getAtomicNum() != Element::S) return false;

  // bonded to exactly three oxygens
  return get_bond_count(mol, atom, Element::O) == 3;
}

inline bool is_S_in_sulfate(const RDKit::RWMol &mol, const RDKit::Atom &atom) {
  if (atom.getAtomicNum() != Element::S) return false;

  // bonded to exactly four oxygens?
  return get_bond_count(mol, atom, Element::O) == 4;
}

inline bool is_P_in_phosphate(const RDKit::RWMol &mol, const RDKit::Atom &atom) {
  if (atom.getAtomicNum() != Element::P) return false;

  unsigned int total_bonds  = get_bond_count(mol, atom);
  unsigned int oxygen_bonds = get_bond_count(mol, atom, 8);

  // all bonds are with oxygen
  return oxygen_bonds == total_bonds;
}

inline bool is_C_in_guanidine(const RDKit::RWMol &mol, const RDKit::Atom &atom) {
  if (atom.getAtomicNum() != Element::C) return false;

  // three bonds and all to nitrogen
  if (get_bond_count(mol, atom) != 3 || get_bond_count(mol, atom, Element::N) != 3) return false;

  // no. of terminal nitrogens (exactly one non-hydrogen bond)
  unsigned int terminal_N_count = 0;
  for (const auto &bond : mol.atomBonds(&atom)) {
    const RDKit::Atom *nbr_atom = bond->getOtherAtom(&atom);

    if (nbr_atom->getAtomicNum() == Element::N) {
      int bonds_to_non_H = get_bond_count(mol, *nbr_atom) - get_bond_count(mol, *nbr_atom, 1);

      if (bonds_to_non_H == 1) {
        terminal_N_count++;
      }
    }
  }

  // two terminal nitrogens
  return terminal_N_count == 2;
}

inline bool is_C_in_acetamidine(const RDKit::RWMol &mol, const RDKit::Atom &atom) {
  if (atom.getAtomicNum() != Element::C) return false;

  // not the most efficient, but reads well
  if (get_bond_count(mol, atom)             != 3 ||
      get_bond_count(mol, atom, Element::N) != 2 ||
      get_bond_count(mol, atom, Element::C) != 1) {
    return false;
  }

  // no. of terminal nitrogens (exactly one non-hydrogen bond)
  unsigned int terminal_N_count = 0;
  for (const auto &bond : mol.atomBonds(&atom)) {
    const RDKit::Atom *nbr = bond->getOtherAtom(&atom);

    if (nbr->getAtomicNum() == Element::N) {
      auto bonds_to_non_H = get_bond_count(mol, *nbr) - get_bond_count(mol, *nbr, 1);

      if (bonds_to_non_H == 1) {
        terminal_N_count++;
      }
    }
  }

  // two terminal nitrogens
  return terminal_N_count == 2;
}

} // namespace lahuta

#endif // LAHUTA_HEURISTICS_HPP
