#include <rdkit/GraphMol/AtomIterators.h>
#include <rdkit/GraphMol/BondIterators.h>

#include "ob/clean_mol.hpp"

namespace lahuta {

namespace {
constexpr std::array<int, 118 + 1> max_bonds{
    0, 1, 0, 1, 2, 4, 4, 4, 2, 1, 0, 1, 2, 6, 6, 6, 6, 1, 0, 1, 2, 6, 6, 6, 6, 8, 6, 6, 6, 6,
    6, 3, 4, 3, 2, 1, 0, 1, 2, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 3, 4, 3, 2, 1, 0, 1, 2, 2, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 3, 4, 3, 2, 1, 0, 1, 2, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
};
} // namespace

// NOTE: Performance is dependent on bond removal (which is the most expensive
// operation). RDKit provides a batch removal utility, which I'm not using.
void clean_bonds(RDKit::RWMol &mol, RDKit::Conformer &conf) {
  RDKit::Bond *maxbond, *bond;
  bool changed;

  for (const auto &atom : mol.atoms()) {
    auto max_bond = max_bonds[atom->getAtomicNum()];
    while (ob_explicit_valence(mol, atom) > max_bond || smallest_bond_angle(mol, conf, atom) < 45.0) {

      maxbond = nullptr;

      // Find first bond
      for (const auto &bond : mol.atomBonds(atom)) {
        maxbond = bond;
        break;
      }

      // No bonds to process
      if (maxbond == nullptr) break;

      // Delete bonds between hydrogens when over max valence
      if (atom->getAtomicNum() == 1) { // Hydrogen
        changed = false;
        for (const auto &bond : mol.atomBonds(atom)) {
          // Neighboring Hydrogen
          if (bond->getOtherAtom(atom)->getAtomicNum() == 1) {
            mol.removeBond(bond->getBeginAtomIdx(), bond->getEndAtomIdx());
            changed = true;
            break;
          }
        }
        if (changed) continue;
      }

      auto max_length = bond_length_sq(conf, maxbond);
      int i = 0;
      for (const auto &bond : mol.atomBonds(atom)) {
        if (i == 0) {
          i++;
          continue;
        }
        auto length = bond_length_sq(conf, bond);
        if (length > max_length) {
          maxbond = bond;
          max_length = length;
        }
      }

      mol.removeBond(maxbond->getBeginAtomIdx(), maxbond->getEndAtomIdx());
    }
  }
}

} // namespace lahuta
