#include <limits>

#include "GraphMol/MonomerInfo.h"
#include "bonds.hpp"
#include "conv.hpp"
#include "nsgrid.hpp"

#include "external/ob/elements.h"

#include "bond_table/bonds.hpp"
#include "bond_table/table.hpp"

using namespace gemmi;

bool connectVdW(Atom p, Atom q, double dist_sq, double covFactor) {
  double rcov1 = covFactor * RDKit::PeriodicTable::getTable()->getRcovalent(
                                 p.element.atomic_number());
  double rcov2 = covFactor * RDKit::PeriodicTable::getTable()->getRcovalent(
                                 q.element.atomic_number());
  if (dist_sq <= (rcov1 + rcov2) * (rcov1 + rcov2)) {
    return true;
  }
  return false;
}

static bool validAdditionalBond(RDKit::Atom *a, RDKit::Atom *b) {
  if (a->getExplicitValence() == 5 && b->getAtomicNum() == 15) {

    // only allow ochedral bonding for F and Cl
    if (a->getAtomicNum() == 9 || a->getAtomicNum() == 17) {
      return true;
    } else {
      return false;
    }
    // other things to check?
    return true;
  }

  if (a->getAtomicNum() == 1 && b->getAtomicNum() == 1) {
    return false;
  }
  if (a->getAtomicNum() == 1 && b->getAtomicNum() == 6) {
    return false;
  }
  if (a->getAtomicNum() == 6 && b->getAtomicNum() == 1) {
    return false;
  }
  return true;
}

bool connectVdW(RDKit::Atom p, RDKit::Atom q, double dist_sq,
                double covFactor) {
  if (dist_sq < 0.16) {
    return false;
  }
  double _rcov1 = covFactor * RDKit::PeriodicTable::getTable()->getRcovalent(
                                  p.getAtomicNum());
  double _rcov2 = covFactor * RDKit::PeriodicTable::getTable()->getRcovalent(
                                  q.getAtomicNum());
  double rcov1 = covFactor * OBElements::GetCovalentRad(p.getAtomicNum());
  double rcov2 = covFactor * OBElements::GetCovalentRad(q.getAtomicNum());
  // std::cout << "Distance: " << dist_sq << std::endl;
  // std::cout << "Rcov1: " << rcov1 << std::endl;
  // std::cout << "Rcov2: " << rcov2 << std::endl;
  // std::cout << "Rcov1*: " << _rcov1 << std::endl;
  // std::cout << "Rcov2*: " << _rcov2 << std::endl;
  if (dist_sq <= (rcov1 + rcov2) * (rcov1 + rcov2)) {
    return true;
  }
  return false;
}

bool connectOBMol(RDKit::Atom *p, RDKit::Atom *q, double dist_sq,
                  double tolerance) {
  if (dist_sq < 0.16) {
    return false;
  }
  // double _rcov1 = tolerance * RDKit::PeriodicTable::getTable()->getRcovalent(
  //                                p.getAtomicNum());
  // double _rcov2 = tolerance * RDKit::PeriodicTable::getTable()->getRcovalent(
  //                                q.getAtomicNum());
  double rcov1 = OBElements::GetCovalentRad(p->getAtomicNum());
  double rcov2 = OBElements::GetCovalentRad(q->getAtomicNum());
  // std::cout << "Distance: " << dist_sq << std::endl;
  // std::cout << "Rcov1: " << rcov1 << std::endl;
  // std::cout << "Rcov2: " << rcov2 << std::endl;
  // std::cout << "Rcov1*: " << _rcov1 << std::endl;
  // std::cout << "Rcov2*: " << _rcov2 << std::endl;
  if (dist_sq <= (rcov1 + rcov2 + tolerance) * (rcov1 + rcov2 + tolerance)) {
    return true;
  }
  return false;
}
bool shouldSkip(RDKit::Atom &atom) {
  auto explicitValence = atom.getExplicitValence();
  if (atom.getExplicitValence() >=
      OBElements::GetMaxBonds(atom.getAtomicNum())) {
    return true;
  }
  if (atom.getAtomicNum() == 7 && atom.getFormalCharge() == 0 &&
      atom.getExplicitValence() >= 3) {
    return true;
  }
  return false;
}

void perceiveBonds(RDKit::RWMol &mol, RDKit::RWMol &newMol,
                   std::vector<int> atomIndices, const NSResults &results,
                   const float covFactor) {

  for (auto i = 0; i < results.getNeighbors().size(); i++) {
    auto res = results.getNeighbors()[i];
    auto dist_sq = results.distances[i];

    auto *a = mol.getAtomWithIdx(res.first);
    auto *b = mol.getAtomWithIdx(res.second);

    // NOTE: first let's handle only within non-protein atoms
    auto aIt = std::find(atomIndices.begin(), atomIndices.end(), a->getIdx());
    auto bIt = std::find(atomIndices.begin(), atomIndices.end(), b->getIdx());
    if (aIt != atomIndices.end() && bIt != atomIndices.end()) {
      if (connectOBMol(a, b, dist_sq, covFactor)) {
        // we need to add the bond to newMol, so we need to get the index of
        // the atoms in newMol that correspond to a and b
        int indexA = std::distance(atomIndices.begin(), aIt);
        int indexB = std::distance(atomIndices.begin(), bIt);
        newMol.addBond(indexA, indexB, RDKit::Bond::BondType::SINGLE);
      }
    }
  }
}

void perceiveBonds(RDKit::RWMol &mol, const NSResults &results,
                   const float covFactor) {

  for (auto i = 0; i < results.getNeighbors().size(); i++) {
    auto res = results.getNeighbors()[i];
    auto dist_sq = results.distances[i];

    auto *a = mol.getAtomWithIdx(res.first);
    auto *b = mol.getAtomWithIdx(res.second);

    // FIXME: Needs to be tested if it is necessary, if so when
    // if (shouldSkip(*a) || shouldSkip(*b)) {
    //   continue;
    // }

    if (connectOBMol(a, b, dist_sq, covFactor)) {
      // FIX: Neighbor indices are guaranteed to be unique, so we don't need
      // this.
      if (mol.getBondBetweenAtoms(a->getIdx(), b->getIdx()) != nullptr) {
        return;
      }
      mol.addBond(a->getIdx(), b->getIdx(), RDKit::Bond::BondType::SINGLE);
    }
  }
}

inline bool is_same_conformer(std::string altlocA, std::string altlocB) {
  return altlocA.empty() || altlocB.empty() || altlocA == altlocB;
}

void logAtomInfoIf(RDKit::Atom *a, RDKit::Atom *b, int idx) {
  if (a->getIdx() == idx || b->getIdx() == idx) {
    std::cout << "Atom 1: " << a->getIdx() << " " << a->getSymbol() << " "
              << a->getMonomerInfo()->getName() << " " << a->getAtomicNum()
              << " " << "Atom 2: " << b->getIdx() << " " << b->getSymbol()
              << " " << b->getMonomerInfo()->getName() << " "
              << b->getAtomicNum() << std::endl;
  }
}

// NOTE: The current approach is to split the molecule into two, one with the
// protein atoms and the other with the non-protein atoms. This split is
// arbitrary, and done only because of predefiend bond orders for protein atoms.
// There are other molecules for which we can define bond orders (water being
// the most important, for example). This would drastically improve performance,
// since the current bottleneck is in SMARTS matching (done for non-protein
// atoms).
// TODO: Implement a struct to hold both RWMol and the non-protein indices
// the non-protein indices should perhaps be put into a set first, to avoid
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

    auto *infoA = static_cast<RDKit::AtomPDBResidueInfo *>(a->getMonomerInfo());
    auto *infoB = static_cast<RDKit::AtomPDBResidueInfo *>(b->getMonomerInfo());

    auto aIsProtein = getToken(infoA->getResidueName());
    auto bIsProtein = getToken(infoB->getResidueName());

    // most likely both atoms are protein atoms, so this check is first
    if (aIsProtein && bIsProtein) {
      if (!is_same_conformer(infoA->getAltLoc(), infoB->getAltLoc())) {
        continue;
      }

      double thresholdB = getElementThreshold(b->getAtomicNum());
      double thresholdA = getElementThreshold(a->getAtomicNum());

      // FIX: Can be further optimized by getting the square directly
      double pairingThreshold = getPairingThreshold(
          a->getAtomicNum(), b->getAtomicNum(), thresholdA, thresholdB);

      if (dist_sq <= pairingThreshold * pairingThreshold) {
        // NOTE: we do not need to check if the bond already exists, since
        // the neighbor indices are guaranteed to be unique
        int order =
            get_intra_bond_order(infoA->getResidueName(), &(infoA->getName()),
                                 &(b->getMonomerInfo()->getName()));
        mol.addBond(a->getIdx(), b->getIdx(), (RDKit::Bond::BondType)order);
      }
      continue;

    } else if (!aIsProtein && !bIsProtein) {

      // handle wanter (TIP3)  special case
      if (infoA->getResidueName() == "TIP3" && infoB->getResidueName() == "TIP3") {
        if (connectOBMol(a, b, dist_sq, 0.45)) {
          mol.addBond(a->getIdx(), b->getIdx(), RDKit::Bond::SINGLE);
        }
        continue;
      }

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
      // if (connectOBMol(a, b, dist_sq, 0.45)) {}
      continue;
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
