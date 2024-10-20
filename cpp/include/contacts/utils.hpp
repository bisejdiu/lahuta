#ifndef LAHUTA_CONTACT_UTILS_HPP
#define LAHUTA_CONTACT_UTILS_HPP

#include <GraphMol/Atom.h>
#include <GraphMol/RWMol.h>

namespace lahuta {

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

} // namespace lahuta

#endif // LAHUTA_CONTACT_UTILS_HPP
