#include "GraphMol/Atom.h"
#include "GraphMol/Bond.h"
#include "GraphMol/MonomerInfo.h"
#include <gemmi/neighbor.hpp>
#include <model.hpp>
#include <unordered_set>
// #include "GraphMol/RWMol.h"
// #include "GraphMol/PeriodicTable.h"
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

void findBondsDeconstructed(Structure &st, Model &model, RDKit::RWMol &mol,
                            double maxRadius, double covFactor) {
  NeighborSearch ns(model, st.cell, maxRadius);
  ns.populate();

  std::vector<int> flags, order, key;

  for (const Chain &chain : model.chains) {
    for (const Residue &res : chain.residues) {
      for (const Atom &atom : res.atoms) {
        // auto indices = model.get_indices(&chain, &res, &atom);
        Position pos = atom.pos;
        char altloc = atom.altloc;

        ns.for_each(
            pos, altloc, maxRadius,
            [&](const NeighborSearch::Mark &m, double dist_sq) {
              const_CRA cra1 = {&chain, &res, &atom};
              const_CRA cra2 = m.to_cra(model);

              if (cra1.atom == cra2.atom) {
                return;
              }

              // skip B-A bonds and only leave A-B bonds
              // if (cra1.atom->serial >= cra2.atom->serial) {
              //   return;
              // }

              if (connectVdW(*cra1.atom, *cra2.atom, dist_sq, covFactor)) {
                // int order = 1;
                // ret.push_back({cra1, cra2, m.image_idx, dist_sq, order});
                // std::cout << "Adding bond between " << atom_str(cra1) << "
                // and " << atom_str(cra2) << std::endl; this is failing for
                // some reason, use try catch to see what's going on Only add if
                // bond already exists
                if (mol.getBondBetweenAtoms(cra1.atom->serial - 1,
                                            cra2.atom->serial - 1) != nullptr) {
                  return;
                }
                mol.addBond((unsigned int)cra1.atom->serial - 1,
                            (unsigned int)cra2.atom->serial - 1,
                            RDKit::Bond::BondType::SINGLE);
              }
            });
      }
    }
  }
};

// std::vector<BondInfo> _findBondsDeconstructed(Structure &st, Model &model,
//                                               double maxRadius) {
//   NeighborSearch ns(model, st.cell, maxRadius);
//   ns.populate();
//
//   std::vector<int> flags, order, key;
//   std::vector<BondInfo> ret;
//
//   initialize_bond_order_table();
//
//   for (const Chain &chain : model.chains) {
//     for (const Residue &res : chain.residues) {
//       for (const gemmi::Atom &atom : res.atoms) {
//         auto indices = model.get_indices(&chain, &res, &atom);
//         Position pos = atom.pos;
//         char altloc = atom.altloc;
//         double thresholdA =
//         getElementThreshold(atom.element.atomic_number());
//
//         ns.for_each(
//             pos, altloc, maxRadius,
//             [&](const NeighborSearch::Mark &m, double dist_sq) {
//
//               const_CRA cra1 = {&chain, &res, &atom};
//               const_CRA cra2 = m.to_cra(model);
//
//               // skip B-A bonds and only leave A-B bonds
//               if (cra1.atom->serial >= cra2.atom->serial) {
//                 return;
//               }
//
//               // ret.push_back({cra1, cra2, m.image_idx, dist_sq});
//
//               const gemmi::Atom *neighborAtom = cra2.atom;
//
//               int neighborElemIdx =
//               getElementIdx(neighborAtom->element.name());
//               // neighborAtom.element.atomic_number();
//               double thresholdB = getElementThreshold(neighborElemIdx);
//               //
//               double dist = std::sqrt(dist_sq);
//               double pairingThreshold =
//                   getPairingThreshold(atom.element.atomic_number(),
//                                       neighborElemIdx, thresholdA,
//                                       thresholdB);
//               //
//               if (dist <= pairingThreshold) {
//                 int order = get_intra_bond_order(res.name, atom.name,
//                                                  neighborAtom->name);
//
//                 ret.push_back({cra1, cra2, m.image_idx, dist_sq, order});
//               }
//             });
//       }
//     }
//   }
//   return ret;
// };
//
//
// //

// inline bool is_same_conformer(std::string altlocA, std::string altlocB) {
//   std::cout << "Comparing " << altlocA << " and " << altlocB << " Should
//   give: " << (altlocA == altlocB) << std::endl; return altlocA == altlocB;
// }
inline bool is_same_conformer(std::string altlocA, std::string altlocB) {
  return altlocA.empty() || altlocB.empty() || altlocA == altlocB;
}

std::vector<int> findBondsDeconstructedRDKit(RDKit::RWMol &mol,
                                             const NSResults &results) {

  std::unordered_set<int> non_protein_indices;

  initialize_bond_order_table();
  // std::cout << "x. size: " << results.getNeighbors().size() << std::endl;
  for (auto i = 0; i < results.getNeighbors().size(); i++) {
    auto res = results.getNeighbors()[i];
    auto dist_sq = results.distances[i];
    auto *a = mol.getAtomWithIdx(res.first);
    auto *b = mol.getAtomWithIdx(res.second);

    auto resiA = a->getMonomerInfo();
    RDKit::AtomPDBResidueInfo *residueA =
        dynamic_cast<RDKit::AtomPDBResidueInfo *>(resiA);
    auto resiB = b->getMonomerInfo();
    RDKit::AtomPDBResidueInfo *residueB =
        dynamic_cast<RDKit::AtomPDBResidueInfo *>(resiB);

    if (!combined_all_names.count(residueA->getResidueName()) &&
        !combined_all_names.count(residueB->getResidueName())) {
      non_protein_indices.insert(a->getIdx());
      non_protein_indices.insert(b->getIdx());
      continue;
    }

    if (!is_same_conformer(residueA->getAltLoc(), residueB->getAltLoc())) {
      continue;
    }

    // if: (1) same chain and (2) same residue and (3) same

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
          get_intra_bond_order(residueA->getResidueName(), residueA->getName(),
                               b->getMonomerInfo()->getName());
      mol.addBond(a->getIdx(), b->getIdx(), (RDKit::Bond::BondType)order);
    }
  }

  std::vector<int> non_protein_indices_vec;
  non_protein_indices_vec.insert(non_protein_indices_vec.end(),
                                 non_protein_indices.begin(),
                                 non_protein_indices.end());

  // sort the indices
  std::sort(non_protein_indices_vec.begin(), non_protein_indices_vec.end());
  return non_protein_indices_vec;
};

RDKit::RWMol lahutaBondAssignment(RDKit::RWMol &mol, const NSResults &results) {

  initialize_bond_order_table();

  std::vector<int> non_protein_indices;
  std::vector<std::pair<int, int>> bonds;

  std::cout << "1\n";
  for (auto i = 0; i < results.getNeighbors().size(); i++) {
    auto res = results.getNeighbors()[i];
    auto dist_sq = results.distances[i];
    auto *a = mol.getAtomWithIdx(res.first);
    auto *b = mol.getAtomWithIdx(res.second);

    auto resiA = a->getMonomerInfo();
    RDKit::AtomPDBResidueInfo *residueA =
        dynamic_cast<RDKit::AtomPDBResidueInfo *>(resiA);
    auto resiB = b->getMonomerInfo();
    RDKit::AtomPDBResidueInfo *residueB =
        dynamic_cast<RDKit::AtomPDBResidueInfo *>(resiB);

    if (!combined_all_names.count(residueA->getResidueName()) &&
        !combined_all_names.count(residueB->getResidueName())) {
      non_protein_indices.push_back(a->getIdx());
      non_protein_indices.push_back(b->getIdx());

      if (connectOBMol(a, b, dist_sq, 0.45)) {
        // bonds.push_back({a->getIdx(), b->getIdx()});
        bonds.emplace_back(a->getIdx(), b->getIdx());
      }
      continue;
    }

    if (!is_same_conformer(residueA->getAltLoc(), residueB->getAltLoc())) {
      continue;
    }

    // if: (1) same chain and (2) same residue and (3) same

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
          get_intra_bond_order(residueA->getResidueName(), residueA->getName(),
                               b->getMonomerInfo()->getName());
      mol.addBond(a->getIdx(), b->getIdx(), (RDKit::Bond::BondType)order);
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
