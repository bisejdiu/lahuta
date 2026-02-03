/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto forward_concat = [](auto&& a, auto&& b, auto&& c) {
 *     return std::string(std::forward<decltype(a)>(a)) +
 *            std::forward<decltype(b)>(b) +
 *            std::forward<decltype(c)>(c);
 *   };
 *   return forward_concat("besian", "sejdiu", "@gmail.com");
 * }();
 *
 */

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>

#include "analysis/sasa/methods/bitmask.hpp"
#include "analysis/sasa/methods/bitmask_simd.hpp"
#include "utils/math_constants.hpp"

namespace lahuta::analysis {

using lahuta::Pi;

template <typename CoordAccessor>
double BitmaskMethod::compute_impl(const ComputeContext<CoordAccessor> &ctx) const {
  const double r_i = ctx.radius();
  if (r_i <= 0.0) return 0.0;

  if (ctx.use_simd) {
    if (lut_.n_points == 64) {
      return simd::compute_bitmask_64_batched(lut_, ctx);
    }
    if (lut_.n_points == 128) {
      return simd::compute_bitmask_128_batched(lut_, ctx);
    }
    if (lut_.n_points == 256) {
      return simd::compute_bitmask_256_batched(lut_, ctx);
    }
  }

  const double r_i_sq     = ctx.radius_sq();
  const std::size_t start = ctx.neighbor_start();
  const std::size_t end   = ctx.neighbor_end();

  // no neighbors means fully exposed
  if (start == end) return 4.0 * Pi * r_i_sq;

  const double inv_two_r = 1.0 / (2.0 * r_i);
  const double cx        = ctx.coords.x(ctx.atom_index);
  const double cy        = ctx.coords.y(ctx.atom_index);
  const double cz        = ctx.coords.z(ctx.atom_index);

  // Initialize visibility bitmask: all points visible
  std::array<std::uint64_t, 4> visibility{};
  for (std::size_t w = 0; w < lut_.words; ++w) {
    visibility[w] = ~std::uint64_t{0};
  }
  if (lut_.words > 0) {
    visibility[lut_.words - 1] = lut_.last_word_mask;
  }

  // Process each neighbor
  for (std::size_t idx = start; idx < end; ++idx) {
    const auto &nbr       = ctx.neighbors.items[idx];
    const std::size_t j   = nbr.index;
    const double vx       = cx - ctx.coords.x(j);
    const double vy       = cy - ctx.coords.y(j);
    const double vz       = cz - ctx.coords.z(j);
    const double v_mag_sq = vx * vx + vy * vy + vz * vz;

    // degen case
    if (v_mag_sq == 0.0) throw std::invalid_argument("SASA invalid geometry: coincident atom centers.");

    const double dist          = std::sqrt(v_mag_sq);
    const double limit         = (nbr.threshold_sq - v_mag_sq - r_i_sq) * inv_two_r;
    const double cos_threshold = limit / dist;

    if (cos_threshold <= -1.0) continue; // no occlusion possible

    // Full occlusion
    if (cos_threshold >= 1.0) {
      for (std::size_t w = 0; w < lut_.words; ++w) {
        visibility[w] = 0;
      }
      break;
    }

    // Look up occlusion mask
    const double inv_dist = 1.0 / dist;
    const int dir_idx   = dir_bin_from_unit(vx * inv_dist, vy * inv_dist, vz * inv_dist, lut_.dir_resolution);
    const int angle_idx = angle_bin_from_cos(cos_threshold, lut_.angle_bins);
    const std::size_t mask_offset = (static_cast<std::size_t>(dir_idx) * lut_.angle_bins + angle_idx);
    const std::uint64_t *mask     = lut_.masks.data() + mask_offset * lut_.words;

    // Apply
    std::uint64_t combined = 0;
    for (std::size_t w = 0; w < lut_.words; ++w) {
      visibility[w] &= ~mask[w];
      combined      |= visibility[w];
    }

    if (combined == 0) break;
  }

  // Count accessible points
  std::size_t accessible = 0;
  for (std::size_t w = 0; w < lut_.words; ++w) {
    accessible += popcount64(visibility[w]);
  }

  const double area_per_point = 4.0 * Pi * r_i_sq / static_cast<double>(lut_.n_points);
  return area_per_point * static_cast<double>(accessible);
}

double BitmaskMethod::compute(const ComputeContext<CoordAccessor> &ctx) const { return compute_impl(ctx); }

template double BitmaskMethod::compute_impl(const ComputeContext<CoordAccessor> &) const;

} // namespace lahuta::analysis
