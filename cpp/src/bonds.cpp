#include "bonds.hpp"
#include "bonds/bonds.hpp"
#include "bonds/table.hpp"
#include "convert.hpp"
#include <rdkit/GraphMol/PeriodicTable.h>
#include <vector>

namespace lahuta {

// FIX: see if this fixes the issue with the periodic table
static const RDKit::PeriodicTable *tbl = RDKit::PeriodicTable::getTable();

// NOTE: this funciton performs two distinct tasks: (1) assign bond orders using
// a table lookup and (2) use smart pattern matching to assign bond orders
// NOTE: simplify by perhaps returning a struct
BondAssignmentResult assign_bonds(RDKit::RWMol &mol, const NSResults &results) {

  std::vector<int> non_predef_atom_indices;
  non_predef_atom_indices.reserve(mol.getNumAtoms());
  std::vector<std::pair<int, int>> bonds;
  std::vector<bool> seen(mol.getNumAtoms(), false);

  std::vector<float> rcov;
  rcov.resize(mol.getNumAtoms());
  for (const auto atom : mol.atoms()) {
    rcov[atom->getIdx()] = tbl->getRcovalent(atom->getAtomicNum());
  }

  for (auto i = 0; i < results.get_pairs().size(); i++) {
    auto res = results.get_pairs()[i];
    auto dist_sq = results.get_distances()[i];
    auto *a = mol.getAtomWithIdx(res.first);
    auto *b = mol.getAtomWithIdx(res.second);

    auto bonded = getIntraBondOrder(a, b);

    if (bonded) { // true only if both atoms are in the table

      double thresholdB = getElementThreshold(b->getAtomicNum());
      double thresholdA = getElementThreshold(a->getAtomicNum());

      // FIX: Get the square directly
      double pairingThreshold = getPairingThreshold(
          a->getAtomicNum(), b->getAtomicNum(), thresholdA, thresholdB);

      int is_a_h = a->getAtomicNum() == 1;
      int is_b_h = b->getAtomicNum() == 1;

      if (dist_sq <= pairingThreshold * pairingThreshold) {
        if (is_a_h ^ is_b_h) {
          auto non_h_atom = a->getAtomicNum() == 1 ? b : a;
          // It is possible to use the following bitwise operation to get the
          // non-hydrogen atom branchless, but it may not lead to a performance
          // gain.
          /*RDKit::Atom* non_h_atom = a + ((b - a) & -(is_a_h));*/
          non_h_atom->setNumExplicitHs(non_h_atom->getNumExplicitHs() + 1);
        }
        mol.addBond(a->getIdx(), b->getIdx(), bonded.bond_type);
      }
      continue;

    } else if (!bonded.atom1_is_predef && !bonded.atom2_is_predef) {

      // NOTE: I am adding bonds regardless of checking the nature of the atoms.
      // This results in metalic bonds being added and likely bonds to different
      // ions, etc.
      if (!seen[a->getIdx()]) {
        non_predef_atom_indices.push_back(a->getIdx());
        seen[a->getIdx()] = true;
      }

      if (!seen[b->getIdx()]) {
        non_predef_atom_indices.push_back(b->getIdx());
        seen[b->getIdx()] = true;
      }

      if (is_bonded_obmol(a, b, dist_sq, 0.45, rcov)) {
        bonds.emplace_back(a->getIdx(), b->getIdx());
      }
      continue;
    } else {

      // if (is_bonded_obmol(a, b, dist_sq, 0.45, rcov)) {
      //   if (a->getSymbol() == "Fe" || b->getSymbol() == "Fe") {
      //     std::cout << "Fe bond: " << a->getIdx() << " " << b->getIdx() <<
      //     std::endl;
      //   mol.addBond(a->getIdx(), b->getIdx(),
      //   RDKit::Bond::BondType::DATIVEONE);
      //   }
      // }

      // non-handled bonds (e.g. matalic bonds, etc.)
      // auto *infoA = static_cast<RDKit::AtomPDBResidueInfo
      // *>(a->getMonomerInfo()); auto *infoB =
      // static_cast<RDKit::AtomPDBResidueInfo *>(b->getMonomerInfo());
      // std::cout << "Unhandled bond: " << a->getIdx() << " " << b->getIdx() <<
      // " " << a->getSymbol() << " " << b->getSymbol() << " " <<
      // infoA->getResidueName() << " " << infoB->getResidueName() << std::endl;
    }
  }

  if (non_predef_atom_indices.empty()) {
    return {};
  }

  std::vector<int> index_mapping;
  index_mapping.resize(mol.getNumAtoms(), -1);
  for (size_t i = 0; i < non_predef_atom_indices.size(); ++i) {
    index_mapping[non_predef_atom_indices[i]] = i;
  }

  auto new_mol = filter_with_conf(mol, non_predef_atom_indices);

  for (const auto &bond : bonds) {
    int aIx = index_mapping[bond.first];
    int bIx = index_mapping[bond.second];

    if (new_mol.getBondBetweenAtoms(aIx, bIx) == nullptr) {
      new_mol.addBond(aIx, bIx, RDKit::Bond::BondType::SINGLE);
    }
  }

  return {new_mol, non_predef_atom_indices};
};

} // namespace lahuta
