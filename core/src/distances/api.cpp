#include <cmath>
#include <stdexcept>
#include <type_traits>

#include "distances/api.hpp"
#include "distances/contiguous_matrix.hpp"
#include "distopia.h"
#include "kd_index.hpp"
#include "nsgrid.hpp"

// clang-format off
namespace lahuta {

namespace detail {
template <typename T>
void flatten_points(const std::vector<std::vector<T>> &points, std::vector<T> &buffer) {
  const size_t n = points.size();
  buffer.resize(3 * n);
  for (size_t i = 0; i < n; ++i) {
    const auto &p = points[i];
    buffer[3 * i + 0] = p[0];
    buffer[3 * i + 1] = p[1];
    buffer[3 * i + 2] = p[2];
  }
}
} // namespace detail

NSResults DistanceComputation::search(const std::vector<std::vector<double>> &points1, const std::vector<std::vector<double>> &points2, double cutoff) {
  RDGeom::POINT3D_VECT tgt;
  tgt.reserve(points1.size());
  for (const auto &r : points1) {
    tgt.emplace_back(r[0], r[1], r[2]);
  }

  RDGeom::POINT3D_VECT qry;
  qry.reserve(points2.size());
  for (const auto &r : points2) {
    qry.emplace_back(r[0], r[1], r[2]);
  }

  KDTreeIndex kd;
  if (!kd.build(tgt)) {
    throw std::runtime_error("KD index build failed: no points");
  }
  return kd.radius_search(qry, cutoff);
}

NSResults DistanceComputation::search(const std::vector<std::vector<double>> &points, double cutoff) {
  FastNS grid(points);
  if (!grid.build(cutoff)) {
    throw std::runtime_error("Box dimension too small for the given cutoff.");
  }
  return grid.self_search();
}

template <typename T>
T DistanceComputation::distance(const std::vector<T> &p1, const std::vector<T> &p2) {
  static_assert(std::is_same<T, float>::value || std::is_same<T, double>::value, "T must be float or double.");

  if (p1.size() < 3 || p2.size() < 3)
    throw std::invalid_argument("Points must have at least 3 dimensions.");

  T d2 = FastNS::dist_sq(p1.data(), p2.data());
  return std::sqrt(d2);
}

template <typename T>
ContiguousMatrix<T> DistanceComputation::distance(const std::vector<std::vector<T>> &points1, const std::vector<std::vector<T>> &points2) {
  static_assert(std::is_same<T, float>::value || std::is_same<T, double>::value, "T must be float or double.");
  const int na = static_cast<int>(points1.size());
  const int nb = static_cast<int>(points2.size());

  if (na == 0 || nb == 0) return ContiguousMatrix<T>(na, nb);

  // Row dim check
  if (points1[0].size() != 3 || points2[0].size() != 3)
    throw std::invalid_argument("Points must have exactly 3 dimensions (x,y,z).");

  std::vector<T> a;
  detail::flatten_points(points1, a);

  std::vector<T> b;
  detail::flatten_points(points2, b);

  ContiguousMatrix<T> res(na, nb);
  distopia::DistanceArrayNoBox(a.data(), b.data(), na, nb, res.data());
  return res;
}

template <typename T>
ContiguousMatrix<T> DistanceComputation::distance(const std::vector<std::vector<T>> &points) {
  static_assert(std::is_same<T, float>::value || std::is_same<T, double>::value, "T must be float or double.");
  const int n = static_cast<int>(points.size());
  ContiguousMatrix<T> result(n, n);

  if (n <= 1) return result;

  std::vector<T> flattened;
  detail::flatten_points(points, flattened);

  const size_t tri_count = static_cast<size_t>(n) * static_cast<size_t>(n - 1) / 2;
  std::vector<T> upper_tri(tri_count);
  distopia::SelfDistanceArrayNoBox(flattened.data(), n, upper_tri.data());

  size_t offset = 0;
  for (int i = 0; i < n; ++i) {
    result(i, i) = static_cast<T>(0);
    for (int j = i + 1; j < n; ++j) {
      const T value = upper_tri[offset++];
      result(i, j) = value;
      result(j, i) = value;
    }
  }

  return result;
}

// Explicit template instantiations
template float DistanceComputation::distance<float>(const std::vector<float> &, const std::vector<float> &);
template double DistanceComputation::distance<double>(const std::vector<double> &, const std::vector<double> &);

template ContiguousMatrix<float> DistanceComputation::distance<float>(const std::vector<std::vector<float>> &, const std::vector<std::vector<float>> &);
template ContiguousMatrix<double> DistanceComputation::distance<double>(const std::vector<std::vector<double>> &, const std::vector<std::vector<double>> &);

template ContiguousMatrix<float> DistanceComputation::distance<float>(const std::vector<std::vector<float>> &);
template ContiguousMatrix<double> DistanceComputation::distance<double>(const std::vector<std::vector<double>> &);

template void detail::flatten_points<float>(const std::vector<std::vector<float>> &, std::vector<float> &);
template void detail::flatten_points<double>(const std::vector<std::vector<double>> &, std::vector<double> &);

} // namespace lahuta
