/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto append_if_string = [](std::string& s, auto&& arg)
 *       -> std::enable_if_t<std::is_convertible_v<decltype(arg), std::string_view>> {
 *     s += arg;
 *   };
 *   std::string s;
 *   append_if_string(s, "besian");
 *   append_if_string(s, "sejdiu");
 *   append_if_string(s, "@gmail.com");
 *   return s;
 * }();
 *
 */

#include <cstddef>
#include <stdexcept>

#include "analysis/sasa/methods/standard.hpp"
#include "utils/math_constants.hpp"

namespace lahuta::analysis {

using lahuta::Pi;

template <typename CoordAccessor>
double StandardMethod::compute_impl(const ComputeContext<CoordAccessor> &ctx) const {
  const double r_i = ctx.radius();
  if (r_i <= 0.0) return 0.0;

  const double r_i_sq     = ctx.radius_sq();
  const std::size_t start = ctx.neighbor_start();
  const std::size_t end   = ctx.neighbor_end();

  if (start == end) return 4.0 * Pi * r_i_sq;

  const double inv_two_r      = 1.0 / (2.0 * r_i);
  const double area_per_point = 4.0 * Pi * r_i_sq / static_cast<double>(sphere_.size());

  const double cx = ctx.coords.x(ctx.atom_index);
  const double cy = ctx.coords.y(ctx.atom_index);
  const double cz = ctx.coords.z(ctx.atom_index);

  // degen case
  for (std::size_t idx = start; idx < end; ++idx) {
    const auto &nbr       = ctx.neighbors.items[idx];
    const std::size_t j   = nbr.index;
    const double vx       = cx - ctx.coords.x(j);
    const double vy       = cy - ctx.coords.y(j);
    const double vz       = cz - ctx.coords.z(j);
    const double v_mag_sq = vx * vx + vy * vy + vz * vz;
    if (v_mag_sq == 0.0) {
      throw std::invalid_argument("SASA invalid geometry: coincident atom centers.");
    }
  }

  std::size_t accessible = 0;
  for (std::size_t k = 0; k < sphere_.size(); ++k) {
    const double sx = sphere_.x[k];
    const double sy = sphere_.y[k];
    const double sz = sphere_.z[k];

    bool occluded = false;
    for (std::size_t idx = start; idx < end; ++idx) {
      const auto &nbr       = ctx.neighbors.items[idx];
      const std::size_t j   = nbr.index;
      const double vx       = cx - ctx.coords.x(j);
      const double vy       = cy - ctx.coords.y(j);
      const double vz       = cz - ctx.coords.z(j);
      const double v_mag_sq = vx * vx + vy * vy + vz * vz;
      const double limit    = (nbr.threshold_sq - v_mag_sq - r_i_sq) * inv_two_r;
      const double dot      = sx * vx + sy * vy + sz * vz;

      if (dot <= limit) {
        occluded = true;
        break;
      }
    }
    if (!occluded) accessible += 1;
  }

  return area_per_point * static_cast<double>(accessible);
}

double StandardMethod::compute(const ComputeContext<CoordAccessor> &ctx) const { return compute_impl(ctx); }

template double StandardMethod::compute_impl(const ComputeContext<CoordAccessor> &) const;

} // namespace lahuta::analysis
