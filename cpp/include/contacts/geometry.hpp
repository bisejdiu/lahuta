/*
    This code adapts functionality from the Mol* (mol-star) molecular visualization package, originally in
   TypeScript, for use in C++ with RDKit in the Lahuta project.

    Mol* is licensed under the MIT License (https://github.com/molstar/molstar), while this adaptation is
   released under the GNU General Public License (GPL).
*/
#ifndef LAHUTA_GEOMETRY_HPP
#define LAHUTA_GEOMETRY_HPP

#include <GraphMol/Atom.h>
#include <GraphMol/RWMol.h>

namespace lahuta {

// Calculate the angles x-a1-a2
// for all x where x is a heavy atom (not H) bonded to atom_a
std::pair<std::vector<double>, std::vector<double>> calculate_angle(
    const RDKit::RWMol &mol, const RDKit::Atom &atom_a, const RDKit::Atom &atom_b,
    bool ignore_hydrogens = true);

double compute_plane_angle(
    const RDGeom::Point3D &atom_a_pos, const RDGeom::Point3D &neighbor1_pos,
    const RDGeom::Point3D &neighbor2_pos, const RDGeom::Point3D &atom_b_pos);

// Function to calculate the angle out of the plane defined by two bonded
// neighbors of atom_a
double compute_plane_angle(const RDKit::RWMol &mol, const RDKit::Atom &atom_a, const RDKit::Atom &atom_b);

} // namespace lahuta

#endif // LAHUTA_GEOMETRY_HPP
