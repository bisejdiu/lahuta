#include <GraphMol/RDKitBase.h> 
#include "ob/clean_mol.hpp"
#include "bonds/table.hpp"

// NOTE: Performance is dependent on bond removal (which is the most expensive
// operation). RDKit provides a batch removal utility, but it is not used here,
// since the function is not a noticeable bottleneck of the application.
void clean_bonds(RDKit::RWMol &mol, RDKit::Conformer &conf) {
  RDKit::Bond *maxbond, *bond;
  bool changed;

  for (auto atomIt = mol.beginAtoms(); atomIt != mol.endAtoms(); ++atomIt) {
    auto atom = *atomIt;
    while (ob_explicit_valence(mol, atom) > max_bonds[atom->getAtomicNum()] ||
           smallest_bond_angle(mol, conf, atom) < 45.0) {

      maxbond = nullptr;

      // Loop through bonds to find the initial maxbond
      for (auto bondIt = mol.getAtomBonds(atom); bondIt.first != bondIt.second;
           ++bondIt.first) {
        bond = mol[*bondIt.first];
        maxbond = bond;
        break;
      }

      if (maxbond == nullptr) {
        break; // No bonds to process
      }

      // Delete bonds between hydrogens when over max valence
      if (atom->getAtomicNum() == 1) { // Hydrogen
        changed = false;
        for (auto bondIt = mol.getAtomBonds(atom);
             bondIt.first != bondIt.second; ++bondIt.first) {
          bond = mol[*bondIt.first];
          // Neighboring Hydrogen
          if (bond->getOtherAtom(atom)->getAtomicNum() == 1) {
            mol.removeBond(bond->getBeginAtomIdx(), bond->getEndAtomIdx());
            changed = true;
            break;
          }
        }
        if (changed) {
          continue; // Reevaluate
        }
      }

      auto maxlength = bond_length_sq(conf, maxbond);
      int i = 0;
      for (auto bondIt = mol.getAtomBonds(atom); bondIt.first != bondIt.second;
           ++bondIt.first) {
        if (i == 0) {
          i++;
          continue;
        }
        bond = mol[*bondIt.first];
        auto length = bond_length_sq(conf, bond);
        if (length > maxlength) {
          maxbond = bond;
          maxlength = length;
        }
      }

      mol.removeBond(maxbond->getBeginAtomIdx(), maxbond->getEndAtomIdx());
    }
  }
}
