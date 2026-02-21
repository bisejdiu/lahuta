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

#include <algorithm>
#include <cmath>
#include <random>
#include <tuple>
#include <vector>

#include <gtest/gtest.h>

#include "distances/api.hpp"
#include "test_utils/fp_test_utils.hpp"

// clang-format off
namespace lahuta {
namespace {

using Point   = std::vector<double>;
using Points  = std::vector<Point>;
using Triplet = std::tuple<int, int, double>;

Points generate_random_points(std::size_t count, std::mt19937 &rng, double min = -500.0, double max = 500.0) {
  using namespace lahuta::test_utils;
  Points points(count, Point(3, 0.0));
  for (auto &p : points) {
    for (double &coord : p) {
      coord = uniform(rng, min, max);
    }
  }
  return points;
}

double brute_distance(const Point &a, const Point &b) {
  double dx = a[0] - b[0];
  double dy = a[1] - b[1];
  double dz = a[2] - b[2];
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}

std::vector<std::vector<double>> brute_distance_matrix(const Points &a, const Points &b) {
  std::vector<std::vector<double>> matrix(a.size(), std::vector<double>(b.size(), 0.0));
  for (std::size_t i = 0; i < a.size(); ++i) {
    for (std::size_t j = 0; j < b.size(); ++j) {
      matrix[i][j] = brute_distance(a[i], b[j]);
    }
  }
  return matrix;
}

std::vector<Triplet> brute_self_neighbors(const Points &points, double cutoff) {
  std::vector<Triplet> neighbors;
  const double cutoff_sq = cutoff * cutoff;
  for (std::size_t i = 0; i < points.size(); ++i) {
    for (std::size_t j = i + 1; j < points.size(); ++j) {
      double dx = points[i][0] - points[j][0];
      double dy = points[i][1] - points[j][1];
      double dz = points[i][2] - points[j][2];
      double d2 = dx * dx + dy * dy + dz * dz;
      if (d2 <= cutoff_sq) {
        neighbors.emplace_back(static_cast<int>(i), static_cast<int>(j), d2);
      }
    }
  }
  std::sort(neighbors.begin(), neighbors.end());
  return neighbors;
}

std::vector<Triplet> brute_cross_neighbors(const Points &reference, const Points &search, double cutoff) {
  std::vector<Triplet> neighbors;
  const double cutoff_sq = cutoff * cutoff;
  for (std::size_t i = 0; i < search.size(); ++i) {
    for (std::size_t j = 0; j < reference.size(); ++j) {
      double dx = search[i][0] - reference[j][0];
      double dy = search[i][1] - reference[j][1];
      double dz = search[i][2] - reference[j][2];
      double d2 = dx * dx + dy * dy + dz * dz;
      if (d2 <= cutoff_sq) {
        neighbors.emplace_back(static_cast<int>(i), static_cast<int>(j), d2);
      }
    }
  }
  std::sort(neighbors.begin(), neighbors.end());
  return neighbors;
}

std::vector<Triplet> canonicalize_self_results(const NSResults &results) {
  std::vector<Triplet> canonical;
  canonical.reserve(results.size());
  const auto &pairs = results.get_pairs();
  const auto &dists = results.get_distances();
  for (std::size_t i = 0; i < pairs.size(); ++i) {
    int first = std::min(pairs[i].first, pairs[i].second);
    int second = std::max(pairs[i].first, pairs[i].second);
    canonical.emplace_back(first, second, static_cast<double>(dists[i]));
  }
  std::sort(canonical.begin(), canonical.end());
  return canonical;
}

std::vector<Triplet> canonicalize_cross_results(const NSResults &results) {
  std::vector<Triplet> canonical;
  canonical.reserve(results.size());
  const auto &pairs = results.get_pairs();
  const auto &dists = results.get_distances();
  for (std::size_t i = 0; i < pairs.size(); ++i) {
    canonical.emplace_back(pairs[i].first, pairs[i].second, static_cast<double>(dists[i]));
  }
  std::sort(canonical.begin(), canonical.end());
  return canonical;
}

} // namespace

TEST(DistanceComputationTest, PairwiseDistanceMatchesBruteForceLargeRandom) {
  std::mt19937 rng(1337);
  const std::size_t count_a = 640;
  const std::size_t count_b = 512;
  auto points_a = generate_random_points(count_a, rng);
  auto points_b = generate_random_points(count_b, rng);

  auto expected = brute_distance_matrix(points_a, points_b);
  auto computed = DistanceComputation::distance(points_a, points_b);

  ASSERT_EQ(computed.rows(), static_cast<int>(points_a.size()));
  ASSERT_EQ(computed.cols(), static_cast<int>(points_b.size()));
  using namespace lahuta::test_utils;
  for (std::size_t i = 0; i < points_a.size(); ++i) {
    for (std::size_t j = 0; j < points_b.size(); ++j) {
      float actual = static_cast<float>(computed(i, j));
      float expect = static_cast<float>(expected[i][j]);
      EXPECT_FLOAT_ULP_EQ(actual, expect, /*max_ulps=*/64, /*abs_guard=*/1e-6f, i, j);
    }
  }
}

TEST(DistanceComputationTest, SelfDistanceMatrixMatchesBruteForce) {
  std::mt19937 rng(2024);
  const std::size_t count = 384;
  auto points = generate_random_points(count, rng, -1000.0, 1000.0);
  auto expected = brute_distance_matrix(points, points);
  auto computed = DistanceComputation::distance(points);

  ASSERT_EQ(computed.rows(), static_cast<int>(points.size()));
  ASSERT_EQ(computed.cols(), static_cast<int>(points.size()));
  using namespace lahuta::test_utils;
  for (std::size_t i = 0; i < points.size(); ++i) {
    for (std::size_t j = 0; j < points.size(); ++j) {
      float actual = static_cast<float>(computed(i, j));
      float expect = static_cast<float>(expected[i][j]);
      EXPECT_FLOAT_ULP_EQ(actual, expect, /*max_ulps=*/64, /*abs_guard=*/1e-6f, i, j);
    }
  }
}

TEST(DistanceComputationTest, SelfDistanceMatrixIsSymmetricAndZeroDiagonal) {
  // Redundant since we have the brute-force oracle, but made explicit for clarity.
  std::vector<std::vector<double>> points{
      {0.0, 0.0, 0.0},
      {1.0, 2.0, 3.0},
      {-4.5, 0.25, 10.0},
      {7.7, -8.8, 9.9},
  };

  auto M = DistanceComputation::distance(points);
  const int n = static_cast<int>(points.size());
  ASSERT_EQ(M.rows(), n);
  ASSERT_EQ(M.cols(), n);

  for (int i = 0; i < n; ++i) {
    EXPECT_DOUBLE_EQ(M(i, i), 0.0) << "Diagonal entry not zero at (" << i << "," << i << ")";
    for (int j = i + 1; j < n; ++j) {
      EXPECT_DOUBLE_EQ(M(i, j), M(j, i)) << "Asymmetry at (" << i << "," << j << ")";
    }
  }
}

TEST(DistanceComputationTest, DistanceBetweenIdenticalPointsIsZero) {
  Point p1{42.0, -17.5, 3.14159};
  Point p2{42.0, -17.5, 3.14159};
  EXPECT_DOUBLE_EQ(DistanceComputation::distance(p1, p2), 0.0);
}

TEST(DistanceComputationTest, DistanceThrowsWhenDimensionIsTooSmall) {
  std::vector<double> bad_point{1.0, 2.0};
  std::vector<double> good_point{1.0, 2.0, 3.0};
  EXPECT_THROW({ DistanceComputation::distance(bad_point, good_point); }, std::invalid_argument);
}

TEST(DistanceComputationTest, DistanceMatrixHandlesEmptyInputs) {
  Points empty;
  Points non_empty = {Point{0.0, 0.0, 0.0}};
  auto result = DistanceComputation::distance(empty, non_empty);
  EXPECT_EQ(result.rows(), 0);
  EXPECT_EQ(result.cols(), 1);

  auto symmetric = DistanceComputation::distance(empty, empty);
  EXPECT_EQ(symmetric.rows(), 0);
  EXPECT_EQ(symmetric.cols(), 0);
}

TEST(DistanceComputationTest, SearchSelfMatchesBruteForceWithDuplicates) {
  Points points = {
      {0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0},
      {1.0, 0.0, 0.0},
      {1.0, 1.0, 1.0},
      {2.0, 2.0, 2.0},
  };

  // Add a random tail to test larger inputs and degeneracies.
  std::mt19937 rng(99);
  auto extra = generate_random_points(128, rng, -5.0, 5.0);
  points.insert(points.end(), extra.begin(), extra.end());

  const double cutoff = 1.5;
  auto expected = brute_self_neighbors(points, cutoff);
  auto results = DistanceComputation::search(points, cutoff);
  auto canonical_results = canonicalize_self_results(results);

  ASSERT_EQ(canonical_results.size(), expected.size());
  for (std::size_t i = 0; i < expected.size(); ++i) {
    EXPECT_EQ  (std::get<0>(canonical_results[i]), std::get<0>(expected[i]));
    EXPECT_EQ  (std::get<1>(canonical_results[i]), std::get<1>(expected[i]));
    EXPECT_NEAR(std::get<2>(canonical_results[i]), std::get<2>(expected[i]), 1e-5);
  }
}

TEST(DistanceComputationTest, SearchCrossMatchesBruteForceOnLargeRandomData) {
  std::mt19937 rng(4242);
  const std::size_t reference_count = 420;
  const std::size_t search_count = 375;

  auto reference_points = generate_random_points(reference_count, rng, -250.0, 250.0);
  auto search_points    = generate_random_points(search_count,    rng, -250.0, 250.0);

  const double cutoff = 120.0;
  auto expected = brute_cross_neighbors(reference_points, search_points, cutoff);
  auto results  = DistanceComputation::search(reference_points, search_points, cutoff);
  auto canonical_results = canonicalize_cross_results(results);

  ASSERT_EQ(canonical_results.size(), expected.size());
  for (std::size_t i = 0; i < expected.size(); ++i) {
    EXPECT_EQ(std::get<0>(canonical_results[i]), std::get<0>(expected[i]));
    EXPECT_EQ(std::get<1>(canonical_results[i]), std::get<1>(expected[i]));
    const double expected_distance_sq = std::get<2>(expected[i]);
    const double tolerance = 1e-3 + 1e-5 * std::abs(expected_distance_sq);
    EXPECT_NEAR(std::get<2>(canonical_results[i]), expected_distance_sq, tolerance);
  }
}

TEST(DistanceComputationTest, SearchHandlesDegenerateDimension) {
  Points points;
  for (int i = 0; i < 5; ++i) {
    points.push_back(Point{static_cast<double>(i), 0.0, 0.0});
  }
  const double cutoff = 1.0;
  auto expected = brute_self_neighbors(points, cutoff);
  auto results  = DistanceComputation::search(points, cutoff);
  auto canonical_results = canonicalize_self_results(results);

  ASSERT_EQ(canonical_results.size(), expected.size());
  for (std::size_t i = 0; i < expected.size(); ++i) {
    EXPECT_EQ  (std::get<0>(canonical_results[i]), std::get<0>(expected[i]));
    EXPECT_EQ  (std::get<1>(canonical_results[i]), std::get<1>(expected[i]));
    EXPECT_NEAR(std::get<2>(canonical_results[i]), std::get<2>(expected[i]), 1e-6);
  }
}

} // namespace lahuta
