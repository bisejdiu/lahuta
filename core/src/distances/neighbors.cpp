/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto count_and_concat = [](auto... args) {
 *     static_assert(sizeof...(args) == 3);
 *     return (std::string{} + ... + std::string(args));
 *   };
 *   return count_and_concat("besian", "sejdiu", "@gmail.com");
 * }();
 *
 */

#include <algorithm>
#include <stdexcept>
#include <vector>

#include <distopia.h>

#include "distances/convert.hpp"
#include "neighbors.hpp"
#include "spatial/fastns.hpp"
#include "spatial/kd_index.hpp"
#include "spatial/nsresults_tls.hpp"

// clang-format off
namespace lahuta::dist {

constexpr int BlockCols = 4096;

namespace {
RDGeom::POINT3D_VECT promote_points(const RDGeom::POINT3D_VECT_F &coords) {
  RDGeom::POINT3D_VECT out;
  out.reserve(coords.size());
  for (const auto &p : coords) {
    out.emplace_back(static_cast<double>(p.x), static_cast<double>(p.y), static_cast<double>(p.z));
  }
  return out;
}
}  // namespace

NSResults neighbors_within_radius_self(const RDGeom::POINT3D_VECT &coords, const NeighborSearchOptions &options) {
  if (coords.empty()) return {};

  FastNS grid(coords);
  if (!grid.build(options.cutoff, options.brute_force_fallback)) {
    if (options.brute_force_fallback) {
      return brute_force_radius_self_streamed(coords, options.cutoff);
    }
    throw std::runtime_error("FastNS grid build failed and brute-force fallback is disabled.");
  }

  return grid.self_search();
}

NSResults neighbors_within_radius_cross(const RDGeom::POINT3D_VECT &queries, const RDGeom::POINT3D_VECT &targets, const NeighborSearchOptions &options) {
  if (queries.empty() || targets.empty()) return {};

  KDTreeIndex index;
  if (!index.build(targets)) {
    auto fallback = brute_force_radius_cross_streamed(queries, targets, options.cutoff);
    return fallback;
  }

  auto results = index.radius_search(queries, options.cutoff);
  return results;
}

NSResults neighbors_within_radius_cross_fastns(const RDGeom::POINT3D_VECT &queries, const RDGeom::POINT3D_VECT &targets, const NeighborSearchOptions &options) {
  if (queries.empty() || targets.empty()) return {};

  FastNS grid(targets);
  if (!grid.build(options.cutoff, options.brute_force_fallback)) {
    if (options.brute_force_fallback) {
      return brute_force_radius_cross_streamed(queries, targets, options.cutoff);
    }
    throw std::runtime_error("FastNS grid build failed and brute-force fallback is disabled.");
  }

  auto results = grid.search(queries);
  return results;
}

NSResults neighbors_within_radius_self(const RDGeom::POINT3D_VECT_F &coords, const NeighborSearchOptions &options) {
  auto promoted = promote_points(coords);
  return neighbors_within_radius_self(promoted, options);
}

NSResults neighbors_within_radius_cross(const RDGeom::POINT3D_VECT_F &queries, const RDGeom::POINT3D_VECT_F &targets, const NeighborSearchOptions &options) {
  auto q = promote_points(queries);
  auto t = promote_points(targets);
  return neighbors_within_radius_cross(q, t, options);
}

NSResults neighbors_within_radius_cross_fastns(const RDGeom::POINT3D_VECT_F &queries, const RDGeom::POINT3D_VECT_F &targets, const NeighborSearchOptions &options) {
  auto q = promote_points(queries);
  auto t = promote_points(targets);
  return neighbors_within_radius_cross_fastns(q, t, options);
}

NSResults brute_force_radius_self_streamed(const RDGeom::POINT3D_VECT &coords, double cutoff) {
  if (coords.empty()) return {};
  const float cutoff_sq = static_cast<float>(cutoff * cutoff);
  auto interleaved = to_interleaved_xyz(coords);

  TlsResultsScope scope;
  NSResults &results = scope.results();
  results.reserve(coords.size());

  std::vector<double> buffer;
  buffer.reserve(BlockCols);

  for (int i = 0; i < static_cast<int>(coords.size()); ++i) {
    const double *ai = &interleaved[3 * i];
    int j = i + 1;
    while (j < static_cast<int>(coords.size())) {
      const int block = std::min(BlockCols, static_cast<int>(coords.size()) - j);
      buffer.resize(block);
      distopia::DistanceArrayNoBox(ai, &interleaved[3 * j], 1, block, buffer.data());
      for (int t = 0; t < block; ++t) {
        const double d = buffer[t];
        const float d2 = static_cast<float>(d * d);
        if (d2 <= cutoff_sq) {
          results.add_neighbors(i, j + t, d2);
        }
      }
      j += block;
    }
  }
  return results;
}

NSResults brute_force_radius_self_streamed(const RDGeom::POINT3D_VECT_F &coords, double cutoff) {
  auto promoted = promote_points(coords);
  return brute_force_radius_self_streamed(promoted, cutoff);
}

NSResults brute_force_radius_cross_streamed(const RDGeom::POINT3D_VECT_F &queries, const RDGeom::POINT3D_VECT_F &targets, double cutoff) {
  auto q = promote_points(queries);
  auto t = promote_points(targets);
  return brute_force_radius_cross_streamed(q, t, cutoff);
}

NSResults brute_force_radius_cross_streamed(const RDGeom::POINT3D_VECT &queries, const RDGeom::POINT3D_VECT &targets, double cutoff) {
  if (queries.empty() || targets.empty()) return {};
  const float cutoff_sq = static_cast<float>(cutoff * cutoff);
  auto query_xyz  = to_interleaved_xyz(queries);
  auto target_xyz = to_interleaved_xyz(targets);

  TlsResultsScope scope;
  NSResults &results = scope.results();
  results.reserve(queries.size());

  std::vector<double> buffer;
  buffer.reserve(BlockCols);

  for (int i = 0; i < static_cast<int>(queries.size()); ++i) {
    const double *qi = &query_xyz[3 * i];
    int j = 0;
    while (j < static_cast<int>(targets.size())) {
      const int block = std::min(BlockCols, static_cast<int>(targets.size()) - j);
      buffer.resize(block);
      distopia::DistanceArrayNoBox(qi, &target_xyz[3 * j], 1, block, buffer.data());
      for (int t = 0; t < block; ++t) {
        const double d = buffer[t];
        const float d2 = static_cast<float>(d * d);
        if (d2 <= cutoff_sq) {
          results.add_neighbors(i, j + t, d2);
        }
      }
      j += block;
    }
  }
  return results;
}

} // namespace lahuta::dist
