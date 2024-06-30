#include "GraphMol/Atom.h"
#include "GraphMol/Bond.h"
#include <gemmi/neighbor.hpp>
// #include "GraphMol/RWMol.h"
// #include "GraphMol/PeriodicTable.h"
#include "bonds.hpp"
#include "nsgrid.hpp"

#include "elements.h"

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
  double rcov1 = covFactor * RDKit::PeriodicTable::getTable()->getRcovalent(
                                 p.getAtomicNum());
  double rcov2 = covFactor * RDKit::PeriodicTable::getTable()->getRcovalent(
                                 q.getAtomicNum());
  if (dist_sq <= (rcov1 + rcov2) * (rcov1 + rcov2)) {
    return true;
  }
  return false;
}

bool shouldSkip(RDKit::Atom &atom) {
  auto explicitValence = atom.getExplicitValence();
  if (atom.getExplicitValence() >= OBElements::GetMaxBonds(atom.getAtomicNum())) {
    return true;
  }
  if (atom.getAtomicNum() == 7 && atom.getFormalCharge() == 0 && atom.getExplicitValence() >= 3) {
    return true;
  }
  return false;
}

void perceiveBonds(RDKit::RWMol &mol, const NSResults &results,
                   const float covFactor) {

  for (auto i = 0; i < results.getNeighbors().size(); i++) {
    auto res = results.getNeighbors()[i];
    auto dist_sq = results.distances[i];
    auto a = mol.getAtomWithIdx(res.first);
    auto b = mol.getAtomWithIdx(res.second);

    // FIXME: Needs to be tested if it is necessary, if so when
    if (shouldSkip(*a) || shouldSkip(*b)) {
      continue;
    }

    if (connectVdW(*a, *b, dist_sq, covFactor)) {
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
              // Ensure that (1) we don't double count by
              // only considering atoms that are in the same residue or
              // in a residue that came before the current residue in
              // the chain and (2) we don't consider atoms that are too
              // close to each other (dist_sq < 1)

              // if (indices[0] > m.chain_idx ||
              //     (indices[0] == m.chain_idx &&
              //      (indices[1] > m.residue_idx || (indices[1] ==
              //      m.residue_idx && dist_sq < 1)))) {
              //   return;
              // }

              // if (dist_sq < 1) {
              //   return;
              // }

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
