/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::array<std::string_view, 4> parts{"besian", "", "sejdiu", "@gmail.com"};
 *   std::vector<std::string_view> valid, empty;
 *   std::partition_copy(parts.begin(), parts.end(), std::back_inserter(valid), std::back_inserter(empty),
 *     [](std::string_view s) { return !s.empty(); });
 *   std::string s; for (auto p : valid) s += p; return s;
 * }();
 *
 */

#ifndef LAHUTA_CONTACTS_HBOND_GEO_VALIDITY_HPP
#define LAHUTA_CONTACTS_HBOND_GEO_VALIDITY_HPP

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>

#include "chemistry/geometry.hpp"
#include "chemistry/neighbors.hpp"
#include "chemistry/utils.hpp"
#include "contacts/molstar/params.hpp"
#include "residues/definitions.hpp"

// clang-format off
namespace lahuta::molstar {

namespace detail {

inline double fast_inv_sqrt(double x) {
  double xhalf = 0.5 * x;
  long long i;
  std::memcpy(&i, &x, sizeof(i));
  i = 0x5fe6ec85e7de30daLL - (i >> 1);
  std::memcpy(&x, &i, sizeof(x));
  x = x * (1.5 - xhalf * x * x);
  return x;
}

inline double safe_inverse_length(const RDGeom::Point3D &vec) {
  const double len_sq = vec.lengthSq();
  if (len_sq <= std::numeric_limits<double>::min()) return 0.0;
  return fast_inv_sqrt(len_sq);
}

inline double normalized_cosine(
  const RDGeom::Point3D &reference,
  double reference_inv_len,
  const RDGeom::Point3D &other
) {
  if (reference_inv_len == 0.0) return 1.0;

  const double len_sq = other.lengthSq();
  if (len_sq <= std::numeric_limits<double>::min()) return 1.0;

  const double inv_len_other = fast_inv_sqrt(len_sq);
  const double cos_val = reference.dotProduct(other) * (reference_inv_len * inv_len_other);
  return std::clamp(cos_val, -1.0, 1.0);
}

} // namespace detail

inline bool is_cluster_sulfur_with_only_metal_neighbors(const RDKit::ROMol &mol, const RDKit::Atom &atom) {
  if (atom.getAtomicNum() != Element::S) return false;

  const auto *res_info = static_cast<const RDKit::AtomPDBResidueInfo *>(atom.getMonomerInfo());
  if (!res_info) return false;
  if (definitions::is_polymer(res_info->getResidueName())) return false;

  bool has_heavy_neighbor = false;
  for (const auto *neighbor : mol.atomNeighbors(&atom)) {
    if (!neighbor || neighbor->getAtomicNum() == 1) continue;
    has_heavy_neighbor = true;
    if (!is_metal_atom(*neighbor)) return false;
  }

  return has_heavy_neighbor;
}

inline bool
are_geometrically_viable(const RDKit::RWMol &mol, const RDKit::Conformer &conf,
                         const RDKit::Atom &donor, const RDKit::Atom &acceptor, const HBondParams &opts) {

  const auto donor_idx          = donor.getIdx();
  const auto acceptor_idx       = acceptor.getIdx();
  const auto &donor_pos         = conf.getAtomPos(donor_idx);
  const auto &acceptor_pos      = conf.getAtomPos(acceptor_idx);
  const auto donor_to_acceptor  = acceptor_pos - donor_pos;
  const double donor_inv_len = detail::safe_inverse_length(donor_to_acceptor);

  const double ideal_don_angle = chemistry::get_atom_geometry_angle(donor.getHybridization());
  const double donor_lower_angle  = std::max(0.0, ideal_don_angle - opts.max_don_angle_dev);
  const double donor_upper_angle  = std::min(M_PI, ideal_don_angle + opts.max_don_angle_dev);
  const double cos_donor_hydrogen = std::cos(opts.max_don_angle_dev);
  const double cos_donor_lower    = std::cos(donor_lower_angle);
  const double cos_donor_upper    = std::cos(donor_upper_angle);
  const bool donor_ignore_metal_neighbors = is_cluster_sulfur_with_only_metal_neighbors(mol, donor);

  bool has_hydrogen_angles  = false;
  bool hydrogen_angle_valid = false;

  for (const auto *neighbor : mol.atomNeighbors(&donor)) {
    if (!neighbor) continue;

    const auto neighbor_idx = neighbor->getIdx();
    if (neighbor_idx == acceptor_idx) continue;

    const auto &neighbor_pos = conf.getAtomPos(neighbor_idx);
    const auto neighbor_vec = neighbor_pos - donor_pos;
    const double cos_angle = detail::normalized_cosine(donor_to_acceptor, donor_inv_len, neighbor_vec);

    if (neighbor->getAtomicNum() == 1) {
      if (opts.ignore_hydrogens) continue;
      has_hydrogen_angles = true;
      if (cos_angle > cos_donor_hydrogen) hydrogen_angle_valid = true;
      continue;
    }

    if (donor_ignore_metal_neighbors && is_metal_atom(*neighbor)) continue;

    if (cos_angle < cos_donor_upper) return false;
    if (donor_lower_angle > 0.0 && cos_angle > cos_donor_lower) return false;
  }

  if (!opts.ignore_hydrogens && has_hydrogen_angles && !hydrogen_angle_valid) return false;

  if (donor.getHybridization() == HybridizationType::SP2) {
    auto out_of_plane = chemistry::compute_plane_angle(mol, conf, donor, acceptor);
    if (out_of_plane.value_or(-1.0) > opts.max_don_out_of_plane_angle) return false;
  }

  const auto &donor_atom_for_acc = (!opts.ignore_hydrogens && has_hydrogen_angles)
      ? chemistry::find_closest_hydrogen_atom(mol, conf, donor, acceptor).value_or(std::cref(donor)).get()
      : donor;

  const auto donor_for_acc_idx = donor_atom_for_acc.getIdx();
  const auto &donor_for_acc_pos = conf.getAtomPos(donor_for_acc_idx);
  const auto acceptor_to_donor = donor_for_acc_pos - acceptor_pos;
  const double acceptor_inv_len = detail::safe_inverse_length(acceptor_to_donor);

  const double ideal_acc_angle = chemistry::get_atom_geometry_angle(acceptor.getHybridization());
  const double min_acc = ideal_acc_angle - opts.max_acc_angle_dev;
  const double cos_min_acc = min_acc <= 0.0 ? 1.0 : std::cos(min_acc);

  for (const auto *neighbor : mol.atomNeighbors(&acceptor)) {
    if (!neighbor) continue;

    const auto neighbor_idx = neighbor->getIdx();
    if (neighbor_idx == donor_for_acc_idx) continue;

    if (neighbor->getAtomicNum() == 1 && opts.ignore_hydrogens) continue;

    const auto &neighbor_pos = conf.getAtomPos(neighbor_idx);
    const auto neighbor_vec = neighbor_pos - acceptor_pos;
    const double cos_angle = detail::normalized_cosine(acceptor_to_donor, acceptor_inv_len, neighbor_vec);

    if (min_acc > 0.0 && cos_angle > cos_min_acc) return false;
  }

  if (acceptor.getHybridization() == RDKit::Atom::HybridizationType::SP2) {
    auto out_of_plane = chemistry::compute_plane_angle(mol, conf, acceptor, donor_atom_for_acc);
    if (out_of_plane.value_or(-1.0) > opts.max_acc_out_of_plane_angle) return false;
  }

  return true;
}

} // namespace lahuta::molstar

#endif // LAHUTA_CONTACTS_HBOND_GEO_VALIDITY_HPP
