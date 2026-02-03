/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
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

#ifndef LAHUTA_MODEL_SSBONDS_HPP
#define LAHUTA_MODEL_SSBONDS_HPP

#include <rdkit/Geometry/point.h>

#include "bonds/rules/table.hpp"
#include "distances/api.hpp"

// clang-format off
namespace lahuta {

inline std::vector<std::pair<int, int>> find_disulfide_bonds(const std::vector<int> &indices, const RDGeom::POINT3D_VECT &coords) {
  std::vector<std::pair<int, int>> result;

  if (indices.size() < 2) return result;

  const double threshold_distance = get_pair_threshold(16, 16);

  size_t n = indices.size();
  size_t num_pairs = n * (n - 1) / 2;

  if (num_pairs < 1000) {
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = i + 1; j < n; ++j) {
        const RDGeom::Point3D &pos1 = coords[indices[i]];
        const RDGeom::Point3D &pos2 = coords[indices[j]];
        double dist = (pos1 - pos2).length();
        if (dist <= threshold_distance) {
          result.emplace_back(indices[i], indices[j]);
        }
      }
    }
  } else {
    // this will almost never get hit
    std::vector<std::vector<double>> sulfur_coords;
    sulfur_coords.reserve(n);
    for (auto sidx : indices) {
      const RDGeom::Point3D &pos = coords[sidx];
      sulfur_coords.push_back({pos.x, pos.y, pos.z});
    }

    auto distances = DistanceComputation::distance(sulfur_coords);
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = i + 1; j < n; ++j) {
        if (distances(i, j) <= threshold_distance) {
          result.emplace_back(indices[i], indices[j]);
        }
      }
    }
  }

  return result;
}

} // namespace lahuta

#endif // LAHUTA_MODEL_SSBONDS_HPP
