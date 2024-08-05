#include <GraphMol/PeriodicTable.h>
#include <GraphMol/RWMol.h>
#include "nsgrid.hpp"

inline bool is_bonded_obmol(RDKit::Atom *p, RDKit::Atom *q, double dist_sq,
                            double tolerance, std::vector<float> &rcov) {
  if (dist_sq < 0.16) {
    return false;
  }

  float rcov_sum = rcov[p->getIdx()] + rcov[q->getIdx()];

  if (dist_sq <= (rcov_sum + tolerance) * (rcov_sum + tolerance)) {
    return true;
  }
  return false;
}

inline bool is_bonded_vdw(RDKit::Atom *p, RDKit::Atom *q, double dist_sq,
                          std::vector<float> &rcov) {
  if (dist_sq < 0.16) {
    return false;
  }

  float rcov_sum = rcov[p->getIdx()] + rcov[q->getIdx()];

  if (dist_sq <= rcov_sum * rcov_sum) {
    return true;
  }
  return false;
}

RDKit::RWMol assign_bonds(RDKit::RWMol &mol, const NSResults &results,
                          std::vector<int> &atom_indices);
