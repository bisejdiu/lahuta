#include "bonds.hpp"
#include "bonds/bonds.hpp"
#include "bonds/table.hpp"
#include "convert.hpp"

auto PeriodicTable = RDKit::PeriodicTable::getTable();

// NOTE: this funciton performs two distinct tasks: (1) assign bond orders using
// a table lookup and (2) use smart pattern matching to assign bond orders
// NOTE: simplify by perhaps returning a struct
RDKit::RWMol assign_bonds(RDKit::RWMol &mol, const NSResults &results,
                          std::vector<int> &non_predef_atom_indices) {

  std::vector<std::pair<int, int>> bonds;
  std::vector<bool> seen(mol.getNumAtoms(), false);

  std::vector<float> rcov;
  rcov.resize(mol.getNumAtoms());
  for (const auto atom : mol.atoms()) {
    rcov[atom->getIdx()] = PeriodicTable->getRcovalent(atom->getAtomicNum());
  }

  for (auto i = 0; i < results.getNeighbors().size(); i++) {
    auto res = results.getNeighbors()[i];
    auto dist_sq = results.distances[i];
    auto *a = mol.getAtomWithIdx(res.first);
    auto *b = mol.getAtomWithIdx(res.second);

    auto order = getIntraBondOrder(a, b);

    if (order) { // true only if both atoms are in the table

      double thresholdB = getElementThreshold(b->getAtomicNum());
      double thresholdA = getElementThreshold(a->getAtomicNum());

      // FIX: Get the square directly
      double pairingThreshold = getPairingThreshold(
          a->getAtomicNum(), b->getAtomicNum(), thresholdA, thresholdB);

      if (dist_sq <= pairingThreshold * pairingThreshold) {
        mol.addBond(a->getIdx(), b->getIdx(),
                    (RDKit::Bond::BondType)(int)order);
      }
      continue;

    } else if (!order.atom1_in_table && !order.atom2_in_table) {

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
      // non-handled bonds (e.g. matalic bonds, etc.)
      // auto *infoA = static_cast<RDKit::AtomPDBResidueInfo
      // *>(a->getMonomerInfo()); auto *infoB =
      // static_cast<RDKit::AtomPDBResidueInfo *>(b->getMonomerInfo());
      // std::cout << "Unhandled bond: " << a->getIdx() << " " << b->getIdx() <<
      // " " << a->getSymbol() << " " << b->getSymbol() << " " <<
      // infoA->getResidueName() << " " << infoB->getResidueName() << std::endl;
    }
  }

  auto newMol = rdMolFromRDKitMol(mol, non_predef_atom_indices);

  std::vector<int> index_mapping;
  index_mapping.resize(mol.getNumAtoms(), -1);

  for (size_t i = 0; i < non_predef_atom_indices.size(); ++i) {
    index_mapping[non_predef_atom_indices[i]] = static_cast<int>(i);
  }

  // Add bonds to the newMol
  for (const auto &bond : bonds) {
    int aIx = index_mapping[bond.first];
    int bIx = index_mapping[bond.second];

    if (newMol.getBondBetweenAtoms(aIx, bIx) == nullptr) {
      newMol.addBond(aIx, bIx, RDKit::Bond::BondType::SINGLE);
    }
  }

  return newMol;
};
