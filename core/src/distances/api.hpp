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

#ifndef LAHUTA_DISTANCES_API_HPP
#define LAHUTA_DISTANCES_API_HPP

#include <vector>

#include "distances/contiguous_matrix.hpp"
#include "spatial/nsresults.hpp"

// clang-format off
namespace lahuta {

// DistanceComputation
// - Element type must be float or double
// - Each point row must have exactly 3 coordinates: x, y, z
// - Containers must be non-empty when a non-empty result is expected
// - Throws std::invalid_argument on malformed input
class DistanceComputation {
public:
  template <typename T>
  static T distance(const std::vector<T> &p1, const std::vector<T> &p2);

  template <typename T>
  static ContiguousMatrix<T> distance(const std::vector<std::vector<T>> &points1, const std::vector<std::vector<T>> &points2);

  template <typename T>
  static ContiguousMatrix<T> distance(const std::vector<std::vector<T>> &points);

  static NSResults search(const std::vector<std::vector<double>> &points1, const std::vector<std::vector<double>> &points2, double cutoff);

  static NSResults search(const std::vector<std::vector<double>> &points, double cutoff);
};

} // namespace lahuta

#endif // LAHUTA_DISTANCES_API_HPP
