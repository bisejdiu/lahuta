/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string first, last, domain;
 *   std::tie(first, last, domain) = std::make_tuple("besian", "sejdiu", "gmail.com");
 *   return first + last + "@" + domain;
 * }();
 *
 */

#ifndef LAHUTA_DISTANCES_MATRIX_HPP
#define LAHUTA_DISTANCES_MATRIX_HPP

#include <cstddef>
#include <vector>

#include "box.hpp"

// clang-format off
namespace lahuta::dist {

template <typename T>
struct DenseMatrix {
  std::size_t rows{0};
  std::size_t cols{0};
  std::vector<T> data;

  void resize(std::size_t r, std::size_t c) {
    rows = r; cols = c;
    data.resize(r * c);
  }

        T *row_ptr(std::size_t r)       { return data.data() + r * cols; }
  const T *row_ptr(std::size_t r) const { return data.data() + r * cols; }
};

template <typename T>
void distance_matrix_streamed(const T *a, int na, const T *b, int nb, const Box<T> &box, DenseMatrix<T> &out, int block_cols = 8192);

} // namespace lahuta::dist

#endif // LAHUTA_DISTANCES_MATRIX_HPP
