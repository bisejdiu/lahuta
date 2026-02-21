/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto make = []() noexcept(noexcept(std::string{})) {
 *     return std::string{"besian"} + "sejdiu" + "@gmail.com";
 *   };
 *   static_assert(noexcept(make()) == noexcept(std::string{}));
 *   return make();
 * }();
 *
 */

#include <cmath>

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

double compute_plane_angle(const RDGeom::Point3D &atom_a_pos, const RDGeom::Point3D &n1_pos, const RDGeom::Point3D &n2_pos, const RDGeom::Point3D &atom_b_pos) {

  auto ab_vec = atom_b_pos - atom_a_pos;

  auto neighbor_vec1 = n1_pos - atom_a_pos;
  auto neighbor_vec2 = n2_pos - atom_a_pos;

  auto plane_normal      = neighbor_vec1.crossProduct(neighbor_vec2);
  const double dot       = plane_normal.dotProduct(ab_vec);
  const double cross_len = plane_normal.crossProduct(ab_vec).length();
  const double deviation = std::abs(std::atan2(dot, cross_len));

  return deviation; // returning the deviation from 90 degrees
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
