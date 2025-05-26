#include "contacts/groups.hpp"
#include "Geometry/point.h"

namespace lahuta {

std::vector<GroupRec> GroupTypeAnalysis::analyze(const RDKit::RWMol &mol, const Residues &residues) {
  std::vector<GroupRec> groups;
  // groups.reserve(...); // ?

  for (std::size_t i = 0; i < NumBuiltinGroups; ++i) {
    auto group_recs = BuilltInGroupFns[i](mol, residues);
    groups.insert(groups.end(), std::make_move_iterator(group_recs.begin()), std::make_move_iterator(group_recs.end()));
  }

  if (!groups.empty() && !groups.front().atoms.empty()) {
    auto &conf = mol.getConformer();

    for (auto &group : groups) {
      RDGeom::Point3D centroid{0.0, 0.0, 0.0};
      for (auto atom_idx : group.atoms) {
          centroid += conf.getAtomPos(atom_idx.get().getIdx());
      }
      centroid /= group.atoms.size();
      group.center = centroid;
    }
  }

  return groups;
}

} // namespace lahuta
