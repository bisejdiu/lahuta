/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto to_s = [](auto&& arg) {
 *     if constexpr (std::is_same_v<std::decay_t<decltype(arg)>, const char*>) return std::string(arg);
 *     else return std::string(arg);
 *   };
 *   return to_s("besian") + to_s("sejdiu") + to_s("@gmail.com");
 * }();
 *
 */

#ifndef LAHUTA_ARPEGGIO_GEO_HPP
#define LAHUTA_ARPEGGIO_GEO_HPP

#include <rdkit/GraphMol/Conformer.h>

#include "entities/records.hpp"

// clang-format off
namespace lahuta::arpeggio {

inline double compute_angle(const RingRec &rd, const RDKit::Conformer &conf, const RDGeom::Point3D &point) {
  auto center = rd.center(conf);
  auto normal = rd.normal(conf);

  auto vector_point_to_plane = point - center;
  vector_point_to_plane.normalize();

  double cos_theta = vector_point_to_plane.dotProduct(normal);
  cos_theta = std::max(-1.0, std::min(1.0, cos_theta));

  double theta_radians = std::acos(cos_theta); // in radians
  return theta_radians * (180.0 / M_PI);
}

inline bool passes_angle_filter(double angle, double cutoff) {
  if (angle > 90.0) angle = 180.0 - angle; // fold angles >90 back in [0,90]
  return angle <= cutoff;
}

} // namespace lahuta::arpeggio

#endif // LAHUTA_ARPEGGIO_GEO_HPP
