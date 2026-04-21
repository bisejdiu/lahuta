/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto noop = [](const char*) {};
 *   std::unique_ptr<const char, decltype(noop)> a("besian", noop);
 *   std::unique_ptr<const char, decltype(noop)> b("sejdiu", noop);
 *   std::unique_ptr<const char, decltype(noop)> c("@gmail.com", noop);
 *   return std::string(a.get()) + b.get() + c.get();
 * }();
 *
 */

/*
    This code adapts functionality from the Mol* (mol-star) molecular visualization package,
    originally in TypeScript, for use in C++ with RDKit in the Lahuta project.

    Mol* is licensed under the MIT License (https://github.com/molstar/molstar),
    while this adaptation is released under the GNU General Public License (GPL).
*/

#ifndef LAHUTA_CHEMISTRY_GEOMETRY_HPP
#define LAHUTA_CHEMISTRY_GEOMETRY_HPP

#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/MonomerInfo.h>
#include <rdkit/GraphMol/RWMol.h>

// clang-format off
namespace lahuta::chemistry {

/// Project a vector onto a plane defined by a normal vector
RDGeom::Point3D project_on_plane(
  const RDGeom::Point3D &vector,
  const RDGeom::Point3D &plane_normal
);

/// Compute the in-plane offset between two points
double compute_in_plane_offset(
  const RDGeom::Point3D &pos_a,
  const RDGeom::Point3D &pos_b,
  const RDGeom::Point3D &normal
);

/// Calculate the angle out of the plane defined by two bonded neighbors of atom_a
double compute_plane_angle(
  const RDGeom::Point3D &atom_a_pos,
  const RDGeom::Point3D &n1_pos,
  const RDGeom::Point3D &n2_pos,
  const RDGeom::Point3D &atom_b_pos
);

/// Function to calculate the angle out of the plane defined by two bonded neighbors of atom_a
std::optional<double> compute_plane_angle(
  const RDKit::RWMol &mol,
  const RDKit::Conformer &conf,
  const RDKit::Atom &atom_a,
  const RDKit::Atom &atom_b
);

} // namespace lahuta::chemistry

#endif // LAHUTA_CHEMISTRY_GEOMETRY_HPP
