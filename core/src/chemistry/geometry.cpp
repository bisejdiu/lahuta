#include "chemistry/geometry.hpp"
#include "chemistry/neighbors.hpp"

// clang-format off
namespace lahuta::chemistry {

RDGeom::Point3D project_on_plane(const RDGeom::Point3D &vector, const RDGeom::Point3D &plane_normal) {
  double scalar_projection = vector.dotProduct(plane_normal);
  RDGeom::Point3D projected_vec = vector - (plane_normal * scalar_projection);
  return projected_vec;
}

double compute_in_plane_offset(const RDGeom::Point3D &pos_a, const RDGeom::Point3D &pos_b, const RDGeom::Point3D &normal) {
  const RDGeom::Point3D vec_ab = pos_a - pos_b;
  const RDGeom::Point3D projected_vec = project_on_plane(vec_ab, normal);
  return projected_vec.length();
}

std::pair<std::vector<double>, std::vector<double>>
calculate_angle(const RDKit::RWMol &mol, const RDKit::Atom &atom_a, const RDKit::Atom &atom_b, bool ignore_hydrogens) {

  std::vector<double> angles;
  std::vector<double> angles_h;

  const auto &conf = mol.getConformer();

  auto &atom_a_pos = conf.getAtomPos(atom_a.getIdx());
  auto &atom_b_pos = conf.getAtomPos(atom_b.getIdx());
  auto ab_vec = atom_b_pos - atom_a_pos;

  for (const auto bond : mol.atomBonds(&atom_a)) {
    auto *other_atom = bond->getOtherAtom(&atom_a);

    if (other_atom->getAtomicNum() == 1 && ignore_hydrogens) continue;

    auto other_atom_pos = conf.getAtomPos(other_atom->getIdx());
    auto ax_vec = other_atom_pos - atom_a_pos;

    double angle = ab_vec.angleTo(ax_vec); // angle between the two vectors (atom_b-atom_a-other_atom)

    (other_atom->getAtomicNum() == 1 ? angles_h : angles).push_back(angle);
  }

  // Return both heavy atom and hydrogen angles
  return std::make_pair(angles, angles_h);
}

double compute_plane_angle(const RDGeom::Point3D &atom_a_pos, const RDGeom::Point3D &n1_pos, const RDGeom::Point3D &n2_pos, const RDGeom::Point3D &atom_b_pos) {

  auto ab_vec = atom_b_pos - atom_a_pos;

  auto neighbor_vec1 = n1_pos - atom_a_pos;
  auto neighbor_vec2 = n2_pos - atom_a_pos;

  auto plane_normal = neighbor_vec1.crossProduct(neighbor_vec2);
  double angle = plane_normal.angleTo(ab_vec);

  return std::abs((M_PI / 2) - angle); // returning the deviation from 90 degrees
}

std::optional<double> compute_plane_angle(const RDKit::RWMol &mol, const RDKit::Atom &atom_a, const RDKit::Atom &atom_b) {

  auto idxs = chemistry::get_bonded_neighbor_indices(atom_a, mol, /*ignore_hydrogens=*/true);
  if (!idxs) return std::nullopt;

  const auto &c = mol.getConformer();
  auto [idx1, idx2] = *idxs;
  auto idx_a = atom_a.getIdx();
  auto idx_b = atom_b.getIdx();
  return compute_plane_angle(c.getAtomPos(idx_a), c.getAtomPos(idx1), c.getAtomPos(idx2), c.getAtomPos(idx_b));
}

} // namespace lahuta::chemistry
