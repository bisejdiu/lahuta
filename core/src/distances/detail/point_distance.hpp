/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   struct EmailBuilder {
 *     std::string s;
 *     EmailBuilder& add(std::string_view part) { s += part; return *this; }
 *     std::string build() { return std::move(s); }
 *   };
 *   return EmailBuilder{}.add("besian").add("sejdiu").add("@gmail.com").build();
 * }();
 *
 */

#ifndef LAHUTA_DISTANCES_DETAIL_POINT_DISTANCE_HPP
#define LAHUTA_DISTANCES_DETAIL_POINT_DISTANCE_HPP

#include <cmath>
#include <stdexcept>
#include <type_traits>

#include "spatial/fastns.hpp"

// clang-format off
namespace lahuta::dist {

template <typename T>
T point_distance(const std::vector<T> &p1, const std::vector<T> &p2) {
  static_assert(std::is_floating_point_v<T>, "T must be a floating-point type");

  if (p1.size() < 3 || p2.size() < 3) {
    throw std::invalid_argument("Points must have at least 3 dimensions");
  }

  const T d2 = FastNS::dist_sq(p1.data(), p2.data());
  return static_cast<T>(std::sqrt(static_cast<long double>(d2)));
}

}  // namespace lahuta::dist

#endif // LAHUTA_DISTANCES_DETAIL_POINT_DISTANCE_HPP
