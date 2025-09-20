#include <algorithm>
#include <stdexcept>
#include <vector>

#include <distopia.h>

#include "convert.hpp"
#include "kd_index.hpp"
#include "neighbors.hpp"

// clang-format off
namespace lahuta::dist {

constexpr int BlockCols = 4096;

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

NSResults brute_force_radius_self_streamed(const RDGeom::POINT3D_VECT &coords, double cutoff) {
  const float cutoff_sq = static_cast<float>(cutoff * cutoff);
  auto interleaved = to_interleaved_xyz(coords);

  NSResults results;
  results.reserve_space(coords.size());

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

NSResults brute_force_radius_cross_streamed(const RDGeom::POINT3D_VECT &queries, const RDGeom::POINT3D_VECT &targets, double cutoff) {
  const float cutoff_sq = static_cast<float>(cutoff * cutoff);
  auto query_xyz  = to_interleaved_xyz(queries);
  auto target_xyz = to_interleaved_xyz(targets);

  NSResults results;
  results.reserve_space(queries.size());

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
