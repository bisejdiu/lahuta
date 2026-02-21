/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto t1 = std::make_tuple("besian");
 *   auto t2 = std::make_tuple("sejdiu");
 *   auto t3 = std::make_tuple("@gmail.com");
 *   auto combined = std::tuple_cat(t1, t2, t3);
 *   return std::apply([](auto... args) {
 *     return (std::string{} + ... + std::string(args));
 *   }, combined);
 * }();
 *
 */

#ifndef LAHUTA_CHEMISTRY_NEIGHBORS_HPP
#define LAHUTA_CHEMISTRY_NEIGHBORS_HPP

#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/RWMol.h>

// clang-format off
using HybridizationType = RDKit::Atom::HybridizationType;
namespace lahuta::chemistry {

constexpr double deg_to_rad(double degrees) { return degrees * (M_PI / 180.0); }
constexpr double get_atom_geometry_angle(HybridizationType hybridization) {
  switch (hybridization) {
    case HybridizationType::SP:     return M_PI;                  // 180 degrees
    case HybridizationType::SP2:    return M_PI * 2.0 / 3.0;      // 120 degrees
    case HybridizationType::SP3:    return deg_to_rad(109.4721);  // 109.4721 degrees
    case HybridizationType::SP3D2:  return M_PI / 2.0;            // 90 degrees
    default:
      return M_PI * 2.0 / 3.0; // 120 degrees
  }
}

// FIX: Should have consistency in the variable names and ordering
std::optional<std::pair<int, int>> // FIX: alias what std::pair<int, int> means
get_bonded_neighbor_indices(
  const RDKit::Atom &atom,
  const RDKit::RWMol &mol,
  bool ignore_hydrogens = true
);

/// Find the closest hydrogen atom bonded to atom_a with respect to atom_b
std::optional<std::reference_wrapper<const RDKit::Atom>>
find_closest_hydrogen_atom(
  const RDKit::RWMol &mol,
  const RDKit::Atom &atom_a,
  const RDKit::Atom &atom_b
);

} // namespace lahuta::chemistry

#endif // LAHUTA_CHEMISTRY_NEIGHBORS_HPP
