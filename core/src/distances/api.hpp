#ifndef LAHUTA_DISTANCES_API_HPP
#define LAHUTA_DISTANCES_API_HPP

#include <vector>

#include "distances/contiguous_matrix.hpp"
#include "nsgrid.hpp"

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
