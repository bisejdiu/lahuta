#include "nsgrid.hpp"
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/RWMol.h>

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

struct BondAssignmentResult {
  RDKit::RWMol mol;
  std::vector<int> atom_indices{};
  bool has_unlisted_resnames = false;

  BondAssignmentResult() = default;
  BondAssignmentResult(RDKit::RWMol mol, std::vector<int> atom_indices)
      : mol(std::move(mol)), atom_indices(std::move(atom_indices)) {
    if (this->atom_indices.size() > 0) {
      this->has_unlisted_resnames = true;
    }
  }
};

BondAssignmentResult assign_bonds(RDKit::RWMol &mol, const NSResults &results);
