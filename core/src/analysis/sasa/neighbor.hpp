/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
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

#ifndef LAHUTA_ANALYSIS_SASA_NEIGHBOR_HPP
#define LAHUTA_ANALYSIS_SASA_NEIGHBOR_HPP

#include <cstddef>
#include <cstdint>
#include <vector>

#include <Geometry/point.h>

#include "analysis/sasa/types.hpp"
#include "distances/neighbors.hpp"
#include "utils/span.hpp"

namespace lahuta::analysis {

// A single neighbor atom.
struct NeighborData {
  AtomIndex index     = 0;
  double threshold_sq = 0.0; // squared expanded radius of the neighbor
};

// Compressed neighbor list for all atoms.
// For atom i, neighbors are in items[offsets[i]..offsets[i+1]).
struct NeighborList {
  std::vector<std::size_t> offsets; // size n+1
  std::vector<NeighborData> items;

  [[nodiscard]] std::size_t atom_count() const noexcept {
    return offsets.empty() ? 0 : offsets.size() - 1; //
  }

  [[nodiscard]] std::size_t neighbor_count(AtomIndex atom_index) const noexcept {
    return offsets[atom_index + 1] - offsets[atom_index];
  }
};

// Coordinate accessor.
struct CoordAccessor {
  span<const RDGeom::Point3D> coords;

  [[nodiscard]] std::size_t size() const noexcept { return coords.size(); }
  [[nodiscard]] double x(AtomIndex i) const noexcept { return coords[i].x; }
  [[nodiscard]] double y(AtomIndex i) const noexcept { return coords[i].y; }
  [[nodiscard]] double z(AtomIndex i) const noexcept { return coords[i].z; }
};

namespace detail {

template <typename Accessor>
double dist_sq(const Accessor &coords, AtomIndex i, AtomIndex j) {
  const double dx = coords.x(i) - coords.x(j);
  const double dy = coords.y(i) - coords.y(j);
  const double dz = coords.z(i) - coords.z(j);
  return dx * dx + dy * dy + dz * dz;
}

} // namespace detail

// Build neighbor lists for SASA computation.
//  - coords Coordinate accessor
//  - expanded_radii Atom radii + probe radius
//  - atom_ids Optional atom IDs for skip_same_id filtering
//  - skip_same_id If true, skip pairs with the same atom_id
template <typename CoordAccessor>
[[nodiscard]] NeighborList build_neighbors(const CoordAccessor &coords, span<const double> expanded_radii,
                                           span<const AtomId> atom_ids, bool skip_same_id,
                                           span<const std::uint8_t> valid = {}) {
  const std::size_t n = coords.size();
  NeighborList list;
  list.offsets.assign(n + 1, 0);
  if (n == 0) return list;

  auto is_valid = [&valid](std::size_t i) noexcept { return valid.empty() || valid[i] != 0; };

  double max_r = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    if (!is_valid(i)) continue;
    const double r = expanded_radii[i];
    if (r > max_r) max_r = r;
  }
  if (max_r <= 0.0) return list;

  const double eps    = std::max(1e-6, 1e-6 * max_r);
  const double cutoff = 2.0 * max_r + eps;

  dist::NeighborSearchOptions opts;
  opts.cutoff               = cutoff;
  opts.brute_force_fallback = true;
  opts.sort_output          = false;

  RDGeom::POINT3D_VECT coord_vec;
  coord_vec.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    coord_vec.emplace_back(coords.x(i), coords.y(i), coords.z(i));
  }

  const auto results = dist::neighbors_within_radius_self(coord_vec, opts);
  const auto &pairs  = results.get_pairs();

  std::vector<std::size_t> counts(n, 0);
  for (const auto &pair : pairs) {
    const auto i = static_cast<std::size_t>(pair.first);
    const auto j = static_cast<std::size_t>(pair.second);
    if (i >= n || j >= n) continue;
    if (!is_valid(i) || !is_valid(j)) continue;
    if (skip_same_id && !atom_ids.empty() && atom_ids[i] == atom_ids[j]) {
      continue;
    }
    const double r_sum = expanded_radii[i] + expanded_radii[j];
    if (r_sum <= 0.0) continue;
    const double d2 = detail::dist_sq(coords, i, j);
    if (d2 <= r_sum * r_sum) {
      counts[i] += 1;
      counts[j] += 1;
    }
  }

  for (std::size_t i = 0; i < n; ++i) {
    list.offsets[i + 1] = list.offsets[i] + counts[i];
  }
  list.items.resize(list.offsets.back());

  std::vector<std::size_t> cursor = list.offsets;
  for (const auto &pair : pairs) {
    const auto i = static_cast<std::size_t>(pair.first);
    const auto j = static_cast<std::size_t>(pair.second);
    if (i >= n || j >= n) continue;
    if (!is_valid(i) || !is_valid(j)) continue;
    if (skip_same_id && !atom_ids.empty() && atom_ids[i] == atom_ids[j]) {
      continue;
    }
    const double r_sum = expanded_radii[i] + expanded_radii[j];
    if (r_sum <= 0.0) continue;
    const double d2 = detail::dist_sq(coords, i, j);
    if (d2 > r_sum * r_sum) continue;

    list.items[cursor[i]++] = NeighborData{j, expanded_radii[j] * expanded_radii[j]};
    list.items[cursor[j]++] = NeighborData{i, expanded_radii[i] * expanded_radii[i]};
  }

  return list;
}

} // namespace lahuta::analysis

#endif // LAHUTA_ANALYSIS_SASA_NEIGHBOR_HPP
