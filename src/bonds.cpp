#include <GraphMol/PeriodicTable.h>
#include "bonds.hpp"
#include "convert.hpp"
#include "bonds/bonds.hpp"
#include "bonds/table.hpp"

auto PeriodicTable = RDKit::PeriodicTable::getTable();

bool connectOBMol(RDKit::Atom *p, RDKit::Atom *q, double dist_sq,
                  double tolerance) {
  if (dist_sq < 0.16) {
    return false;
  }
  double rcov1 = PeriodicTable->getRcovalent(p->getAtomicNum());
  double rcov2 = PeriodicTable->getRcovalent(q->getAtomicNum());

  if (dist_sq <= (rcov1 + rcov2 + tolerance) * (rcov1 + rcov2 + tolerance)) {
    return true;
  }
  return false;
}

// TODO: Implement a struct to hold both RWMol and the non-protein indices the
// non-protein indices should perhaps be put into a set first, to avoid
// duplicates.
// Function to collect unique indices during iteration over neighboring pairs
RDKit::RWMol lahutaBondAssignment(RDKit::RWMol &mol, const NSResults &results,
                                  std::vector<int> &non_protein_indices) {

  std::vector<std::pair<int, int>> bonds;
  std::vector<bool> seen(mol.getNumAtoms(), false);

  constexpr int MAX_INDEX = std::numeric_limits<int>::max();
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
        non_protein_indices.push_back(a->getIdx());
        seen[a->getIdx()] = true;
      }

      if (!seen[b->getIdx()]) {
        non_protein_indices.push_back(b->getIdx());
        seen[b->getIdx()] = true;
      }

      if (connectOBMol(a, b, dist_sq, 0.45)) {
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

  auto newMol = rdMolFromRDKitMol(mol, non_protein_indices);

  std::vector<int> index_mapping;
  index_mapping.resize(mol.getNumAtoms(), -1);

  for (size_t i = 0; i < non_protein_indices.size(); ++i) {
    index_mapping[non_protein_indices[i]] = static_cast<int>(i);
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
