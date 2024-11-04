#include "contacts/geometry.hpp"
#include "contacts/hydrogen_bonds.hpp"

namespace lahuta {
namespace geometry {

std::pair<std::vector<double>, std::vector<double>> calculate_angle(
    const RDKit::RWMol &mol, const RDKit::Atom &atom_a, const RDKit::Atom &atom_b, bool ignore_hydrogens) {

  std::vector<double> angles;
  std::vector<double> angles_h;

  const auto &conf = mol.getConformer();

  auto &atom_a_pos = conf.getAtomPos(atom_a.getIdx());
  auto &atom_b_pos = conf.getAtomPos(atom_b.getIdx());
  auto ab_vec = atom_b_pos - atom_a_pos;

  for (const auto bond : mol.atomBonds(&atom_a)) {
    auto *other_atom = bond->getOtherAtom(&atom_a);

    // possibly ignore hydrogens
    if (other_atom->getAtomicNum() == 1 && ignore_hydrogens) {
      continue;
    }

    auto other_atom_pos = conf.getAtomPos(other_atom->getIdx());
    auto ax_vec = other_atom_pos - atom_a_pos;

    // angle between the two vectors (atom_b-atom_a-other_atom)
    double angle = ab_vec.angleTo(ax_vec);

    (other_atom->getAtomicNum() == 1 ? angles_h : angles).push_back(angle);
  }

  // Return both heavy atom and hydrogen angles
  return std::make_pair(angles, angles_h);
}

double compute_plane_angle(
    const RDGeom::Point3D &atom_a_pos, const RDGeom::Point3D &neighbor1_pos,
    const RDGeom::Point3D &neighbor2_pos, const RDGeom::Point3D &atom_b_pos) {

  auto ab_vec = atom_b_pos - atom_a_pos;

  auto neighbor_vec1 = neighbor1_pos - atom_a_pos;
  auto neighbor_vec2 = neighbor2_pos - atom_a_pos;

  // plane normal
  auto plane_normal = neighbor_vec1.crossProduct(neighbor_vec2);

  // angle between plane's normal vector and ab_vec
  double angle = plane_normal.angleTo(ab_vec);

  // returning the deviation from 90 degrees
  return std::abs((M_PI / 2) - angle);
}

double compute_plane_angle(const RDKit::RWMol &mol, const RDKit::Atom &atom_a, const RDKit::Atom &atom_b) {

  const RDKit::Conformer &conf = mol.getConformer();
  auto &atom_a_pos = conf.getAtomPos(atom_a.getIdx());
  auto &atom_b_pos = conf.getAtomPos(atom_b.getIdx());

  // Retrieve up to two neighbor positions
  auto neighbor_positions = get_neighbor_positions(atom_a, conf, mol);

  if (neighbor_positions.size() < 2) {
    return -1.0;
  }

  // Compute the angle using the positions
  return compute_plane_angle(atom_a_pos, *neighbor_positions[0], *neighbor_positions[1], atom_b_pos);
}

} // namespace geometry
} // namespace lahuta
