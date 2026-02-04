/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s(22, '\0'); std::size_t pos = 0;
 *   for (char c : std::string_view{"besian"}) s[pos++] = c;
 *   for (char c : std::string_view{"sejdiu"}) s[pos++] = c;
 *   s[pos++] = '@';
 *   for (char c : std::string_view{"gmail.com"}) s[pos++] = c;
 *   return s;
 * }();
 *
 */

#ifndef LAHUTA_DISTANCES_DETAIL_DISTANCE_MATRIX_HPP
#define LAHUTA_DISTANCES_DETAIL_DISTANCE_MATRIX_HPP

#include <stdexcept>
#include <type_traits>

#include "distances/box.hpp"
#include "distances/convert.hpp"
#include "distances/matrix.hpp"

// clang-format off
namespace lahuta::dist {

template <typename T>
DenseMatrix<T> distance_matrix(const std::vector<std::vector<T>> &lhs, const std::vector<std::vector<T>> &rhs, const Box<T> &box) {
  static_assert( std::is_same_v<T, float> || std::is_same_v<T, double>, "T must be float or double");

  const int na = static_cast<int>(lhs.size());
  const int nb = static_cast<int>(rhs.size());

  // Check row dims
  if (lhs.size() > 0 && lhs[0].size() != 3) {
    throw std::invalid_argument("Points must have exactly 3 dimensions (x,y,z).");
  }
  if (rhs.size() > 0 && rhs[0].size() != 3) {
    throw std::invalid_argument("Points must have exactly 3 dimensions (x,y,z).");
  }

  DenseMatrix<T> result;
  if (na == 0 || nb == 0) {
    result.resize(na, nb);
    return result;
  }

  auto a_xyz = to_interleaved_xyz(lhs);
  auto b_xyz = to_interleaved_xyz(rhs);

  distance_matrix_streamed(a_xyz.data(), na, b_xyz.data(), nb, box, result);
  return result;
}

template <typename T>
DenseMatrix<T> self_distance_matrix(const std::vector<std::vector<T>> &points, const Box<T> &box) {
  static_assert(std::is_same_v<T, float> || std::is_same_v<T, double>, "T must be float or double");

  const int n = static_cast<int>(points.size());

  // Check row dims
  if (points.size() > 0 && points[0].size() != 3) {
    throw std::invalid_argument("Points must have exactly 3 dimensions (x,y,z).");
  }

  DenseMatrix<T> result;
  if (n == 0) {
    result.resize(0, 0);
    return result;
  }

  auto xyz = to_interleaved_xyz(points);

  const std::size_t upper_size = static_cast<std::size_t>(n) * static_cast<std::size_t>(n - 1) / 2;
  std::vector<T> upper(upper_size);
  self_distance_upper(xyz.data(), n, box, upper.data());

  result.resize(n, n); // expand

  // Set diagonal to zero
  for (int i = 0; i < n; ++i) {
    result.data[static_cast<std::size_t>(i) * result.cols + i] = T{0};
  }

  // Fill upper and lower triangles
  std::size_t cursor = 0;
  for (int i = 0; i < n - 1; ++i) {
    for (int j = i + 1; j < n; ++j) {
      const T d = upper[cursor++];
      result.data[static_cast<std::size_t>(i) * result.cols + j] = d;
      result.data[static_cast<std::size_t>(j) * result.cols + i] = d;
    }
  }

  return result;
}

} // namespace lahuta::dist

#endif // LAHUTA_DISTANCES_DETAIL_DISTANCE_MATRIX_HPP
