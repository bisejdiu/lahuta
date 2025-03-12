#ifndef LAHUTA_DISTANCES_HPP
#define LAHUTA_DISTANCES_HPP

#include <cmath>
#include <stdexcept>
#include <vector>

#include "distopia.h"
#include "nsgrid.hpp"

#include <type_traits>

namespace lahuta {

template <typename T> using Vector = std::vector<T>;
template <typename T> using Matrix = std::vector<std::vector<T>>;

template<typename T>
class ContiguousMatrix {
public:
    ContiguousMatrix(int rows, int cols) : rows_(rows), cols_(cols), data_(rows * cols, T()) {}

    T&       operator()(int i, int j)       { return data_[i * cols_ + j]; }
    const T& operator()(int i, int j) const { return data_[i * cols_ + j]; }

    T*       data()       { return data_.data(); }
    const T* data() const { return data_.data(); }

    int rows() const { return rows_; }
    int cols() const { return cols_; }

private:
    int rows_, cols_;
    std::vector<T> data_;
};

class DistanceComputation {
public:
  template <typename T>
  static T distance(const Vector<T> &p1, const Vector<T> &p2) {
    static_assert(std::is_floating_point_v<T>, "T must be a floating-point type.");

    if (p1.size() < 3 || p2.size() < 3) {
      throw std::invalid_argument("Points must have at least 3 dimensions.");
    }
    T d2 = FastNS::dist_sq(p1.data(), p2.data());
    return std::sqrt(d2);
  }


  // FIX: use distopia to speed things up
  template <typename T>
  static Matrix<T> distance_old(const Matrix<T> &points1, const Matrix<T> &points2) {
    size_t n1 = points1.size();
    size_t n2 = points2.size();
    Matrix<T> result(n1, Vector<T>(n2, 0.0));
    for (size_t i = 0; i < n1; ++i) {
      for (size_t j = 0; j < n2; ++j) {
        result[i][j] = distance(points1[i], points2[j]);
      }
    }
    return result;
  }

  template <typename T>
  static ContiguousMatrix<T> distance(const std::vector<std::vector<T>> &points1, const std::vector<std::vector<T>> &points2) {
    const int na = static_cast<int>(points1.size());
    const int nb = static_cast<int>(points2.size());

    if (na == 0 || nb == 0) return ContiguousMatrix<T>(na, nb);

    std::vector<T> a(3 * na);
    for (int i = 0; i < na; ++i) {
      a[3 * i + 0] = points1[i][0];
      a[3 * i + 1] = points1[i][1];
      a[3 * i + 2] = points1[i][2];
    }

    std::vector<T> b(3 * nb);
    for (int i = 0; i < nb; ++i) {
      b[3 * i + 0] = points2[i][0];
      b[3 * i + 1] = points2[i][1];
      b[3 * i + 2] = points2[i][2];
    }

    ContiguousMatrix<T> res(na, nb);
    distopia::DistanceArrayNoBox(a.data(), b.data(), na, nb, res.data());
    return res;
  }

  template <typename T> 
  static ContiguousMatrix<T> distance(const Vector<T> &points) { return distance(points, points); }

  // Compute distances between two sets of points given a cutoff.
  static NSResults search(const Matrix<double> &points1, const Matrix<double> &points2, double cutoff) {
    FastNS grid(points1);
    if (!grid.build(cutoff)) {
      throw std::runtime_error("Box dimension too small for the given cutoff.");
    }
    return grid.search(points2);
  }

  static NSResults search(const Matrix<double> &points, double cutoff) {
    FastNS grid(points);
    if (!grid.build(cutoff)) {
      throw std::runtime_error("Box dimension too small for the given cutoff.");
    }
    return grid.self_search();
  }
};

} // namespace lahuta

#endif // LAHUTA_DISTANCES_HPP
