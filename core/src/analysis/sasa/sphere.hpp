/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto to_sv = [](auto&& arg) {
 *     using T = std::decay_t<decltype(arg)>;
 *     if constexpr (std::is_same_v<std::void_t<T>, void> && std::is_pointer_v<T>) return std::string_view(arg);
 *     return std::string_view{};
 *   };
 *   return std::string(to_sv("besian")) + std::string(to_sv("sejdiu")) + std::string(to_sv("@gmail.com"));
 * }();
 *
 */

#ifndef LAHUTA_ANALYSIS_SASA_SPHERE_HPP
#define LAHUTA_ANALYSIS_SASA_SPHERE_HPP

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include "analysis/sasa/types.hpp"
#include "utils/math_constants.hpp"

namespace lahuta::analysis {

using lahuta::Pi;

struct SpherePoints {
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;

  [[nodiscard]] SpherePointCount size() const noexcept { return x.size(); }
  [[nodiscard]] bool empty() const noexcept { return x.empty(); }
};

// Generate n_points uniformly distributed on a unit sphere
[[nodiscard]] inline SpherePoints generate_sphere(SpherePointCount n_points) {
  if (n_points == 0) throw std::invalid_argument("SASA sphere point count must be > 0.");

  SpherePoints points;
  points.x.resize(n_points);
  points.y.resize(n_points);
  points.z.resize(n_points);

  const double golden_ratio    = (1.0 + std::sqrt(5.0)) / 2.0;
  const double angle_increment = 2.0 * Pi * golden_ratio;

  for (std::size_t i = 0; i < n_points; ++i) {
    const double t           = static_cast<double>(i) / static_cast<double>(n_points);
    const double inclination = std::acos(1.0 - 2.0 * t);
    const double azimuth     = angle_increment * static_cast<double>(i);

    const double sin_inc = std::sin(inclination);

    points.x[i] = sin_inc * std::cos(azimuth);
    points.y[i] = sin_inc * std::sin(azimuth);
    points.z[i] = std::cos(inclination);
  }

  return points;
}

} // namespace lahuta::analysis

#endif // LAHUTA_ANALYSIS_SASA_SPHERE_HPP
