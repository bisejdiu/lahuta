#ifndef LAHUTA_CONTACT_UTILS_HPP
#define LAHUTA_CONTACT_UTILS_HPP

#include "GraphMol/MonomerInfo.h"
#include <GraphMol/Atom.h>
#include <GraphMol/RWMol.h>

namespace lahuta {

constexpr double deg_to_rad(double degrees) { return degrees * (M_PI / 180.0); }

/// Returns the atoms bonded to the given atom with the specified atomic number.
inline std::vector<const RDKit::Atom *>
bonded_atoms(const RDKit::RWMol &mol, const RDKit::Atom *atom, int atomic_number) {
  std::vector<const RDKit::Atom *> bonded_atoms;
  for (const auto &bond : mol.atomBonds(atom)) {
    RDKit::Atom *neighbor_atom = bond->getOtherAtom(atom);
    if (neighbor_atom->getAtomicNum() == atomic_number) {
      bonded_atoms.push_back(neighbor_atom);
    }
  }
  return bonded_atoms;
}

/// Returns the number of bonds for a given atom in the molecule.
inline unsigned int get_bond_count(const RDKit::RWMol &mol, const RDKit::Atom &atom) {
  unsigned int bond_count = 0;
  for (const auto &bond : mol.atomBonds(&atom)) {
    bond_count++;
  }
  return bond_count;
}

/// Returns the number of bonds where the neighboring atom has the specified atomic number.
inline unsigned int
get_bond_count(const RDKit::RWMol &mol, const RDKit::Atom &atom, unsigned int atomic_num) {
  unsigned int bond_count = 0;
  for (const auto &bond : mol.atomBonds(&atom)) {
    const RDKit::Atom *neighbor_atom = bond->getOtherAtom(&atom);
    if (neighbor_atom->getAtomicNum() == atomic_num) {
      bond_count++;
    }
  }
  return bond_count;
}

/// Returns the number of hydrogen atoms bonded to the given atom.
inline int get_h_count(RDKit::ROMol &mol, RDKit::Atom &atom) {
  int hCount = 0;
  for (const auto &bondIt : mol.atomBonds(&atom)) {
    auto other_atom = bondIt->getOtherAtom(&atom);
    if (other_atom->getAtomicNum() == 1) {
      hCount++;
    }
  }
  return hCount;
}

inline bool is_same_residue(const RDKit::RWMol &mol, const RDKit::Atom &atom_a, const RDKit::Atom &atom_b) {
  auto info_a = static_cast<const RDKit::AtomPDBResidueInfo *>(atom_a.getMonomerInfo());
  auto info_b = static_cast<const RDKit::AtomPDBResidueInfo *>(atom_b.getMonomerInfo());

  if (!info_a || !info_b) return false;

  return info_a->getResidueNumber() == info_b->getResidueNumber()
         && info_a->getResidueName() == info_b->getResidueName()
         && info_a->getChainId() == info_b->getChainId();
}

// FIX: need a is_protein check
inline bool are_residueids_close(
    const RDKit::RWMol &mol, const RDKit::Atom &atom_a, const RDKit::Atom &atom_b, int threshold) {
  auto info_a = static_cast<const RDKit::AtomPDBResidueInfo *>(atom_a.getMonomerInfo());
  auto info_b = static_cast<const RDKit::AtomPDBResidueInfo *>(atom_b.getMonomerInfo());

  if (!info_a || !info_b) return false;
  if (info_a->getChainId() != info_b->getChainId()) return false;

  return std::abs(info_a->getResidueNumber() - info_b->getResidueNumber()) <= threshold;
}

} // namespace lahuta

#endif // LAHUTA_CONTACT_UTILS_HPP
