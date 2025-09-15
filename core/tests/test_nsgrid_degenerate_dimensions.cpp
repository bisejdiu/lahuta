#include <gtest/gtest.h>

#include <vector>

#include "nsgrid.hpp"

// clang-format off
namespace lahuta {

// These tests verify that FastNS handles degenerate
// dimensions correctly without stalling or creating massive grids.

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
  float scale_factor = 1.1f;

  FastNS ns(coords, scale_factor);

  bool build_ok = ns.build(cutoff);
  EXPECT_FALSE(build_ok) << "Build should fail for degenerate data";
}

TEST(NSGridDegenerateTest, CompletelyDegenerateData) {
  std::vector<std::vector<double>> coords = {
      {-3.47773801e-107, -3.47773801e-107, -3.47773801e-107},
      {-3.47773801e-107, -3.47773801e-107, -3.47773801e-107},
      {-3.47773801e-107, -3.47773801e-107, -3.47773801e-107},
      {-3.47773801e-107, -3.47773801e-107, -3.47773801e-107},
      {-3.47773801e-107, -3.47773801e-107, -3.47773801e-107}};

  double cutoff = 1.0;
  FastNS ns(coords, 1.1f);

  bool build_ok = ns.build(cutoff);
  EXPECT_FALSE(build_ok) << "Should fail for completely degenerate data";
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
  FastNS ns(coords, 1.1f);

  bool build_ok = ns.build(cutoff);

  if (build_ok) {
    auto results = ns.self_search();

    // Results cannot exceed n*(n-1)/2 pairs
    const size_t n = coords.size();
    const size_t max_pairs = n * (n - 1) / 2;
    EXPECT_LE(results.size(), max_pairs);

    // Distances are squared and must be <= cutoff**2, indices must be in range and not self-pairs
    const float cutoff2 = static_cast<float>(cutoff * cutoff);
    const auto &pairs = results.get_pairs();
    const auto &dists = results.get_distances();
    ASSERT_EQ(pairs.size(), dists.size());
    for (size_t k = 0; k < dists.size(); ++k) {
      EXPECT_LE(dists[k], cutoff2);
      int i = pairs[k].first;
      int j = pairs[k].second;
      EXPECT_NE(i, j);
      EXPECT_GE(i, 0);
      EXPECT_GE(j, 0);
      EXPECT_LT(i, static_cast<int>(n));
      EXPECT_LT(j, static_cast<int>(n));
    }
  }
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
  FastNS ns(coords, 1.1f);

  bool build_ok = ns.build(cutoff);

  if (build_ok) {
    auto results = ns.self_search();

    const size_t n = coords.size();
    const size_t max_pairs = n * (n - 1) / 2;
    EXPECT_LE(results.size(), max_pairs);

    const float cutoff2 = static_cast<float>(cutoff * cutoff);
    const auto &pairs = results.get_pairs();
    const auto &dists = results.get_distances();
    ASSERT_EQ(pairs.size(), dists.size());
    for (size_t k = 0; k < dists.size(); ++k) {
      EXPECT_LE(dists[k], cutoff2);
      int i = pairs[k].first;
      int j = pairs[k].second;
      EXPECT_NE(i, j);
      EXPECT_GE(i, 0);
      EXPECT_GE(j, 0);
      EXPECT_LT(i, static_cast<int>(n));
      EXPECT_LT(j, static_cast<int>(n));
    }
  }
}

} // namespace lahuta
