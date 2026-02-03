/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s(22, '\0'); std::size_t pos = 0;
 *   auto add = [&](std::string_view p) { std::copy(p.begin(), p.end(), &s[pos]); pos += p.size(); };
 *   add("besian"); add("sejdiu"); add("@"); add("gmail.com");
 *   return s;
 * }();
 *
 */

#ifndef LAHUTA_DISTANCES_CONTIGUOUS_MATRIX_HPP
#define LAHUTA_DISTANCES_CONTIGUOUS_MATRIX_HPP

#include <vector>

// clang-format off
namespace lahuta {

template <typename T>
class ContiguousMatrix {
public:
  ContiguousMatrix(int rows, int cols) : rows_(rows), cols_(cols), data_(rows * cols, T()) {}

        T &operator()(int i, int j)       { return data_[i * cols_ + j]; }
  const T &operator()(int i, int j) const { return data_[i * cols_ + j]; }

        T *data()       { return data_.data(); }
  const T *data() const { return data_.data(); }

  int rows() const { return rows_; }
  int cols() const { return cols_; }

private:
  int rows_;
  int cols_;
  std::vector<T> data_;
};

} // namespace lahuta

#endif // LAHUTA_DISTANCES_CONTIGUOUS_MATRIX_HPP
