#include <algorithm>
#include <cmath>
#include <random>
#include <tuple>
#include <vector>

#include <gtest/gtest.h>

#include "distances/neighbors.hpp"
#include "test_utils/fp_test_utils.hpp"

namespace lahuta::dist {

namespace {
using Triplet = std::tuple<int, int, float>; // (i,j,d2)

static RDGeom::POINT3D_VECT make_random_points(std::size_t n, std::mt19937 &rng, double lo, double hi) {
  using namespace lahuta::test_utils;
  RDGeom::POINT3D_VECT pts;
  pts.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    pts.emplace_back(uniform(rng, lo, hi), uniform(rng, lo, hi), uniform(rng, lo, hi));
  }
  return pts;
}

static std::vector<Triplet> brute_self_oracle(const RDGeom::POINT3D_VECT &coords, double cutoff) {
  const float cutoff_sq = static_cast<float>(cutoff * cutoff);
  std::vector<Triplet> out;
  const int n = static_cast<int>(coords.size());
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      const double dx = coords[i].x - coords[j].x;
      const double dy = coords[i].y - coords[j].y;
      const double dz = coords[i].z - coords[j].z;
      const float d2 = static_cast<float>(dx * dx + dy * dy + dz * dz);
      if (d2 <= cutoff_sq) out.emplace_back(i, j, d2);
    }
  }
  std::sort(out.begin(), out.end());
  return out;
}

static std::vector<Triplet> canonicalize(const NSResults &res) {
  std::vector<Triplet> out;
  out.reserve(res.size());
  const auto &pairs = res.get_pairs();
  const auto &dists = res.get_distances();
  for (std::size_t k = 0; k < pairs.size(); ++k) {
    int i = pairs[k].first, j = pairs[k].second;
    if (j < i) std::swap(i, j);
    out.emplace_back(i, j, dists[k]);
  }
  std::sort(out.begin(), out.end());
  return out;
}

} // namespace

TEST(FastNSGridBehavior, RandomCloudMatchesManual) {
  std::mt19937 rng(13371337);
  const auto coords = make_random_points(750, rng, -100.0, 100.0);
  const double cutoff = 7.5; // keeps neighbor count reasonable

  const auto oracle = brute_self_oracle(coords, cutoff);

  NeighborSearchOptions opts;
  opts.cutoff = cutoff;

  auto results = neighbors_within_radius_self(coords, opts);
  auto canon = canonicalize(results);

  ASSERT_EQ(canon.size(), oracle.size());
  for (std::size_t k = 0; k < oracle.size(); ++k) {
    EXPECT_EQ(std::get<0>(canon[k]), std::get<0>(oracle[k]));
    EXPECT_EQ(std::get<1>(canon[k]), std::get<1>(oracle[k]));
    const float expected = std::get<2>(oracle[k]);
    const float actual = std::get<2>(canon[k]);
    const float tol = 2e-4f + 2e-6f * std::fabs(expected);
    EXPECT_NEAR(actual, expected, tol);
  }
}

TEST(FastNSGridBehavior, BorderlineBoxSizeAboveBelowCutoff) {
  // Create a small cube whose size can be adjusted to be just above/below cutoff
  std::mt19937 rng(4242);
  const double cutoff = 1.0;
  using namespace lahuta::test_utils;

  auto make_cube = [&](double span) {
    RDGeom::POINT3D_VECT pts;
    pts.reserve(300);
    for (int i = 0; i < 300; ++i)
      pts.emplace_back(
          uniform(rng, 0.0, 1.0) * span,
          uniform(rng, 0.0, 1.0) * span,
          uniform(rng, 0.0, 1.0) * span);
    return pts;
  };

  const double eps = 1e-4;           // nudge around the threshold
  const double span_below = cutoff - eps;
  const double span_above = cutoff + eps;

  auto below = make_cube(span_below);
  auto above = make_cube(span_above);

  const auto oracle_below = brute_self_oracle(below, cutoff);
  const auto oracle_above = brute_self_oracle(above, cutoff);

  NeighborSearchOptions opts;
  opts.cutoff = cutoff;

  auto res_below = neighbors_within_radius_self(below, opts);
  auto can_below = canonicalize(res_below);
  ASSERT_EQ(can_below.size(), oracle_below.size());
  for (std::size_t k = 0; k < oracle_below.size(); ++k) {
    EXPECT_EQ(std::get<0>(can_below[k]), std::get<0>(oracle_below[k]));
    EXPECT_EQ(std::get<1>(can_below[k]), std::get<1>(oracle_below[k]));
    const float expected = std::get<2>(oracle_below[k]);
    const float actual = std::get<2>(can_below[k]);
    const float tol = 2e-4f + 2e-6f * std::fabs(expected);
    EXPECT_NEAR(actual, expected, tol);
  }

  auto res_above = neighbors_within_radius_self(above, opts);
  auto can_above = canonicalize(res_above);
  ASSERT_EQ(can_above.size(), oracle_above.size());
  for (std::size_t k = 0; k < oracle_above.size(); ++k) {
    EXPECT_EQ(std::get<0>(can_above[k]), std::get<0>(oracle_above[k]));
    EXPECT_EQ(std::get<1>(can_above[k]), std::get<1>(oracle_above[k]));
    const float expected = std::get<2>(oracle_above[k]);
    const float actual = std::get<2>(can_above[k]);
    const float tol = 2e-4f + 2e-6f * std::fabs(expected);
    EXPECT_NEAR(actual, expected, tol);
  }
}

} // namespace lahuta::dist
