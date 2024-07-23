#include "GraphMol/MonomerInfo.h"
#include "bonds.hpp"
#include "conv.hpp"
#include "nsgrid.hpp"

#include "elements.h"

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

RDKit::RWMol lahutaBondAssignment(RDKit::RWMol &mol, const NSResults &results,
                                  std::vector<int> &non_protein_indices) {

  initialize_bond_order_table();

  std::vector<std::pair<int, int>> bonds;

  std::cout << "1\n";
  for (auto i = 0; i < results.getNeighbors().size(); i++) {
    auto res = results.getNeighbors()[i];
    auto dist_sq = results.distances[i];
    auto *a = mol.getAtomWithIdx(res.first);
    auto *b = mol.getAtomWithIdx(res.second);

    auto *infoA = static_cast<RDKit::AtomPDBResidueInfo *>(a->getMonomerInfo());
    auto *infoB = static_cast<RDKit::AtomPDBResidueInfo *>(b->getMonomerInfo());

    bool aIsProtein = combined_all_names.find(infoA->getResidueName()) !=
                      combined_all_names.end();
    bool bIsProtein = combined_all_names.find(infoB->getResidueName()) !=
                      combined_all_names.end();

    // most likely both atoms are protein atoms
    if (aIsProtein && bIsProtein) {
      if (!is_same_conformer(infoA->getAltLoc(), infoB->getAltLoc())) {
        continue;
      }

      double thresholdB = getElementThreshold(b->getAtomicNum());
      double thresholdA = getElementThreshold(a->getAtomicNum());

      // these are likely going to be too many lookups. Needs to be optimized
      double pairingThreshold = getPairingThreshold(
          a->getAtomicNum(), b->getAtomicNum(), thresholdA, thresholdB);

      if (dist_sq <= pairingThreshold * pairingThreshold) {
        // if (mol.getBondBetweenAtoms(a->getIdx(), b->getIdx()) != nullptr) {
        //   continue;
        // }
        int order =
            get_intra_bond_order(infoA->getResidueName(), &(infoA->getName()),
                                 &(b->getMonomerInfo()->getName()));
        mol.addBond(a->getIdx(), b->getIdx(), (RDKit::Bond::BondType)order);
      }
      continue;
    } else if (!aIsProtein && !bIsProtein) {
      non_protein_indices.push_back(a->getIdx());
      non_protein_indices.push_back(b->getIdx());
      if (connectOBMol(a, b, dist_sq, 0.45)) {
        logAtomInfoIf(a, b, 60115);
        bonds.emplace_back(a->getIdx(), b->getIdx());
      }
      continue;
    } else {
      // logAtomInfoIf(a, b, 60115);
      if (connectOBMol(a, b, dist_sq, 0.45)) {
        // double thresholdB = getElementThreshold(b->getAtomicNum());
        // double thresholdA = getElementThreshold(a->getAtomicNum());
        //
        // // these are likely going to be too many lookups. Needs to be
        // optimized double pairingThreshold = getPairingThreshold(
        //     a->getAtomicNum(), b->getAtomicNum(), thresholdA, thresholdB);
        // if (dist_sq <= pairingThreshold * pairingThreshold) {
        // std::cout << "Unhandled potential bond: " << a->getIdx() << " - "
        //           << b->getIdx() << std::endl;
        // logAtomInfoIf(a, b, 60115);
      }
      continue;
    }
  }

  std::sort(non_protein_indices.begin(), non_protein_indices.end());
  non_protein_indices.erase(
      std::unique(non_protein_indices.begin(), non_protein_indices.end()),
      non_protein_indices.end());

  std::cout << "2\n";
  auto newMol = rdMolFromRDKitMol(mol, non_protein_indices);
  std::cout << "3\n";
  // Map old indices to new indices
  std::unordered_map<int, int> old_to_new_index;
  for (size_t i = 0; i < non_protein_indices.size(); ++i) {
    old_to_new_index[non_protein_indices[i]] = i;
  }

  // Add bonds to the newMol
  for (const auto &bond : bonds) {
    auto it_a = old_to_new_index.find(bond.first);
    auto it_b = old_to_new_index.find(bond.second);

    if (it_a != old_to_new_index.end() && it_b != old_to_new_index.end()) {
      int aIx = it_a->second;
      int bIx = it_b->second;

      if (newMol.getBondBetweenAtoms(aIx, bIx) == nullptr) {
        newMol.addBond(aIx, bIx, RDKit::Bond::BondType::SINGLE);
      }
    }
  }
  std::cout << "4\n";

  return newMol;
};
