/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s{"besian"};
 *   s.append("sejdiu").append("@gmail.com");
 *   return s;
 * }();
 *
 */

#include <algorithm>
#include <cmath>
#include <tuple>
#include <vector>

#include <gtest/gtest.h>

#include "spatial/fastns.hpp"
#include "spatial/nsresults.hpp"

// clang-format off
namespace lahuta {

// These tests verify that FastNS handles degenerate
// dimensions correctly without stalling or creating massive grids.
namespace {

using Triplet = std::tuple<int, int, float>;

std::vector<Triplet> brute_force_pairs(const std::vector<std::vector<double>> &coords, double cutoff) {
  const int n = static_cast<int>(coords.size());
  const float cutoff_sq = static_cast<float>(cutoff * cutoff);
  std::vector<Triplet> out;
  out.reserve(static_cast<std::size_t>(n) * static_cast<std::size_t>(n - 1) / 2);
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      const double dx = coords[i][0] - coords[j][0];
      const double dy = coords[i][1] - coords[j][1];
      const double dz = coords[i][2] - coords[j][2];
      const float d2 = static_cast<float>(dx * dx + dy * dy + dz * dz);
      if (d2 <= cutoff_sq) out.emplace_back(i, j, d2);
    }
  }
  std::sort(out.begin(), out.end());
  return out;
}

std::vector<Triplet> canonicalize(const NSResults &results) {
  std::vector<Triplet> out;
  const auto &pairs = results.get_pairs();
  const auto &dists = results.get_distances();
  out.reserve(pairs.size());
  for (std::size_t k = 0; k < pairs.size(); ++k) {
    int i = pairs[k].first;
    int j = pairs[k].second;
    if (j < i) std::swap(i, j);
    out.emplace_back(i, j, dists[k]);
  }
  std::sort(out.begin(), out.end());
  return out;
}

void expect_fastns_matches_brute(const std::vector<std::vector<double>> &coords, double cutoff) {
  FastNS ns(coords);
  ASSERT_TRUE(ns.build(cutoff)) << "FastNS build should succeed";
  auto actual   = canonicalize(ns.self_search());
  auto expected = brute_force_pairs(coords, cutoff);
  ASSERT_EQ(actual.size(), expected.size());
  for (std::size_t k = 0; k < expected.size(); ++k) {
    EXPECT_EQ(std::get<0>(actual[k]), std::get<0>(expected[k])) << k;
    EXPECT_EQ(std::get<1>(actual[k]), std::get<1>(expected[k])) << k;
    const float exp = std::get<2>(expected[k]);
    const float act = std::get<2>(actual[k]);
    const float tol = 1e-3f + 1e-5f * std::fabs(exp);
    EXPECT_NEAR(act, exp, tol) << k;
  }
}

} // namespace

// Test data that originally caused infinite stalls in hypothesis tests
TEST(NSGridDegenerateTest, PropertyBasedStallCase) {
  std::vector<std::vector<double>> coords = {
      {-3.47773801e-107, -3.47773801e-107, -3.47773801e-107},
      {-3.47773801e-107, -3.47773801e-107, -3.47773801e-107},
      {-3.47773801e-107, -3.47773801e-107, -3.47773801e-107},
      { 1.41544969e000,  -3.47773801e-107, -3.47773801e-107},
      {-3.48044136e000,  -3.47773801e-107, -3.47773801e-107},
      {-3.47773801e-107, -3.47773801e-107, -3.47773801e-107},
      {-3.47773801e-107, -3.47773801e-107, -3.47773801e-107},
      {-3.47773801e-107, -3.47773801e-107, -3.47773801e-107},
      {-3.47773801e-107, -3.47773801e-107, -3.47773801e-107},
      {-3.47773801e-107, -3.47773801e-107, -3.47773801e-107},
      {-3.47773801e-107, -3.47773801e-107, -3.47773801e-107},
      {-3.47773801e-107, -3.47773801e-107, -3.47773801e-107},
      {-3.47773801e-107, -3.47773801e-107, -5.00000000e000},
      {-3.47773801e-107, -3.47773801e-107, -3.47773801e-107},
      {-9.07279265e-001, -3.47773801e-107, -3.47773801e-107},
      {-3.47773801e-107, -3.47773801e-107, -3.47773801e-107},
      {-3.47773801e-107, -3.47773801e-107, -3.47773801e-107},
      {-3.47773801e-107,  1.34856373e-040, -3.47773801e-107},
      {-3.47773801e-107, -3.47773801e-107, -3.47773801e-107},
      {-3.47773801e-107, -3.47773801e-107, -3.47773801e-107},
      {-4.38939717e-001, -3.47773801e-107, -3.47773801e-107},
      {-1.50965892e-291, -3.47773801e-107, -3.47773801e-107},
      {-3.47773801e-107, -3.47773801e-107, -3.47773801e-107},
      {-3.47773801e-107, -3.47773801e-107, -3.47773801e-107},
      {-3.47773801e-107, -3.47773801e-107, -3.47773801e-107}};

  double cutoff = 2.4956705026482906;
  expect_fastns_matches_brute(coords, cutoff);
}

TEST(NSGridDegenerateTest, CompletelyDegenerateData) {
  std::vector<std::vector<double>> coords = {
      {-3.47773801e-107, -3.47773801e-107, -3.47773801e-107},
      {-3.47773801e-107, -3.47773801e-107, -3.47773801e-107},
      {-3.47773801e-107, -3.47773801e-107, -3.47773801e-107},
      {-3.47773801e-107, -3.47773801e-107, -3.47773801e-107},
      {-3.47773801e-107, -3.47773801e-107, -3.47773801e-107}};

  double cutoff = 1.0;
  expect_fastns_matches_brute(coords, cutoff);
}

TEST(NSGridDegenerateTest, MixedNormalAndDegenerateDimensions) {
  std::vector<std::vector<double>> coords = {
      {1.41544969e000,   -3.47773801e-107, -3.47773801e-107},
      {-3.48044136e000,  -3.47773801e-107, -3.47773801e-107},
      {-3.47773801e-107, -3.47773801e-107, -5.00000000e000},
      {-9.07279265e-001, -3.47773801e-107, -3.47773801e-107},
      {-3.47773801e-107,  1.34856373e-040, -3.47773801e-107},
      {-4.38939717e-001, -3.47773801e-107, -3.47773801e-107},
      {-1.50965892e-291, -3.47773801e-107, -3.47773801e-107}
  };

  double cutoff = 2.5;
  expect_fastns_matches_brute(coords, cutoff);
}

TEST(NSGridDegenerateTest, NumericalEdgeCases) {
  std::vector<std::vector<double>> coords = {
      {0.0, 0.0, 0.0},
      {1e-100, 1e-100, 1e-100},
      {-1e-100, -1e-100, -1e-100},
      {1e-200, 0.0, 0.0},
      {0.0, 1e-200, 0.0},
      {0.0, 0.0, 1e-200}};

  double cutoff = 1.0;
  expect_fastns_matches_brute(coords, cutoff);
}

} // namespace lahuta
