/*
    This code adapts functionality from the Mol* (mol-star) molecular visualization package,
    originally in TypeScript, for use in C++ with RDKit in the Lahuta project.

    Mol* is licensed under the MIT License (https://github.com/molstar/molstar),
    while this adaptation is released under the GNU General Public License (GPL).
*/

#ifndef LAHUTA_GEOMETRY_HPP
#define LAHUTA_GEOMETRY_HPP

#include "Geometry/point.h"
#include <GraphMol/Atom.h>
#include <GraphMol/RWMol.h>

namespace lahuta {
namespace geometry {

/// Project a vector onto a plane defined by a normal vector
inline RDGeom::Point3D project_on_plane(const RDGeom::Point3D &vector, const RDGeom::Point3D &plane_normal) {
  // subtract component along the normal
  double scalar_projection = vector.dotProduct(plane_normal);
  RDGeom::Point3D projected_vec = vector - (plane_normal * scalar_projection);
  return projected_vec;
}

/// Compute the in-plane offset between two points
inline double compute_in_plane_offset(
    const RDGeom::Point3D &pos_a, const RDGeom::Point3D &pos_b, const RDGeom::Point3D &normal) {
  RDGeom::Point3D vec_ab = pos_a - pos_b;
  RDGeom::Point3D projected_vec = project_on_plane(vec_ab, normal);
  double in_plane_offset = projected_vec.length();
  return in_plane_offset;
}

/// Calculate the angles x-a1-a2 for all x where x is a heavy atom (not H) bonded to atom_a
std::pair<std::vector<double>, std::vector<double>> calculate_angle(
    const RDKit::RWMol &mol, const RDKit::Atom &atom_a, const RDKit::Atom &atom_b,
    bool ignore_hydrogens = true);

/// Function to calculate the angle out of the plane defined by two bonded neighbors of atom_a
double compute_plane_angle(
    const RDGeom::Point3D &atom_a_pos, const RDGeom::Point3D &neighbor1_pos,
    const RDGeom::Point3D &neighbor2_pos, const RDGeom::Point3D &atom_b_pos);

/// Function to calculate the angle out of the plane defined by two bonded neighbors of atom_a
double compute_plane_angle(const RDKit::RWMol &mol, const RDKit::Atom &atom_a, const RDKit::Atom &atom_b);

} // namespace geometry
} // namespace lahuta

#endif // LAHUTA_GEOMETRY_HPP
