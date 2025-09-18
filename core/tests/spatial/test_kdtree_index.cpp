#include <gtest/gtest.h>

#include <algorithm>
#include <random>
#include <tuple>
#include <vector>

#include "distances/neighbors.hpp"
#include "kd_index.hpp"
#include "test_utils/fp_test_utils.hpp"

// clang-format off
namespace lahuta::dist {

namespace {

using Triplet = std::tuple<int, int, float>; // (query_idx, target_idx, d2)

RDGeom::POINT3D_VECT make_random_points(int n, double lo, double hi, std::mt19937 &rng) {
  RDGeom::POINT3D_VECT pts;
  pts.reserve(n);
  for (int i = 0; i < n; ++i) {
    pts.emplace_back(
        test_utils::uniform(rng, lo, hi),
        test_utils::uniform(rng, lo, hi),
        test_utils::uniform(rng, lo, hi));
  }
  return pts;
}

std::vector<Triplet> brute_cross_oracle(const RDGeom::POINT3D_VECT &queries, const RDGeom::POINT3D_VECT &targets, double cutoff) {
  const float cutoff_sq = static_cast<float>(cutoff * cutoff);
  std::vector<Triplet> out;
  const int nq = static_cast<int>(queries.size());
  const int nt = static_cast<int>(targets.size());
  out.reserve(std::min(nq * nt, 200000));
  for (int i = 0; i < nq; ++i) {
    for (int j = 0; j < nt; ++j) {
      const double dx = queries[i].x - targets[j].x;
      const double dy = queries[i].y - targets[j].y;
      const double dz = queries[i].z - targets[j].z;
      const float d2 = static_cast<float>(dx * dx + dy * dy + dz * dz);
      if (d2 <= cutoff_sq) out.emplace_back(i, j, d2);
    }
  }
  std::sort(out.begin(), out.end());
  return out;
}

std::vector<Triplet> canonicalize(const NSResults &res) {
  std::vector<Triplet> out;
  out.reserve(res.size());
  const auto &pairs = res.get_pairs();
  const auto &d2    = res.get_distances();
  for (std::size_t k = 0; k < pairs.size(); ++k) {
    out.emplace_back(pairs[k].first, pairs[k].second, d2[k]);
  }
  std::sort(out.begin(), out.end());
  return out;
}

void expect_match(const std::vector<Triplet> &actual, const std::vector<Triplet> &expected, const char *label) {
  ASSERT_EQ(actual.size(), expected.size()) << label;
  for (std::size_t k = 0; k < expected.size(); ++k) {
    EXPECT_EQ(std::get<0>(actual[k]), std::get<0>(expected[k])) << label;
    EXPECT_EQ(std::get<1>(actual[k]), std::get<1>(expected[k])) << label;
    EXPECT_FLOAT_ULP_EQ(
        std::get<2>(actual[k]),
        std::get<2>(expected[k]),
        /*max_ulps=*/64,
        /*abs_guard=*/1e-6f,
        std::get<0>(expected[k]),
        std::get<1>(expected[k])
    );
  }
}

} // namespace

TEST(KDTreeIndex, CrossMatchesBruteForceRandomPoints) {
  std::mt19937 rng(1337);
  const double cutoff = 2.0;
  auto targets = make_random_points(/*n=*/256, -5.0, 5.0, rng);
  auto queries = make_random_points(/*n=*/180, -5.0, 5.0, rng);

  KDTreeIndex index;
  ASSERT_TRUE(index.build(targets));

  auto kd_results    = canonicalize(index.radius_search(queries, cutoff));
  auto brute_results = brute_cross_oracle(queries, targets, cutoff);
  expect_match(kd_results, brute_results, "random-owned");
}

TEST(KDTreeIndex, BuildViewFloat32MatchesOwnedResults) {
  std::mt19937 rng(424242);
  const int n_targets = 240;
  const double cutoff = 1.5;

  std::vector<float> coords;
  coords.reserve(static_cast<std::size_t>(n_targets) * 3);
  RDGeom::POINT3D_VECT targets;
  targets.reserve(n_targets);
  for (int i = 0; i < n_targets; ++i) {
    const double x = test_utils::uniform(rng, -4.0, 4.0);
    const double y = test_utils::uniform(rng, -4.0, 4.0);
    const double z = test_utils::uniform(rng, -4.0, 4.0);
    coords.push_back(static_cast<float>(x));
    coords.push_back(static_cast<float>(y));
    coords.push_back(static_cast<float>(z));
    targets.emplace_back(x, y, z);
  }

  auto queries = make_random_points(/*n=*/160, -4.0, 4.0, rng);

  KDTreeIndex index;
  ASSERT_TRUE(index.build_view_f32(coords.data(), targets.size(), /*leaf_size=*/24));
  auto kd_results    = canonicalize(index.radius_search(queries, cutoff));
  auto brute_results = brute_cross_oracle(queries, targets, cutoff);
  expect_match(kd_results, brute_results, "view-f32");

  // Rebuild on the same instance with owned storage to ensure state reset works.
  ASSERT_TRUE(index.build(targets));
  auto rebuilt_results = canonicalize(index.radius_search(queries, cutoff));
  expect_match(rebuilt_results, brute_results, "rebuild-owned");
}

TEST(KDTreeIndex, BuildViewFloat64UsesDoubleTree) {
  std::mt19937 rng(2025);
  const int n_targets = 128;
  const double cutoff = 1.2;

  std::vector<double> coords;
  coords.reserve(static_cast<std::size_t>(n_targets) * 3);
  RDGeom::POINT3D_VECT targets;
  targets.reserve(n_targets);
  for (int i = 0; i < n_targets; ++i) {
    const double x = test_utils::uniform(rng, -3.5, 3.5);
    const double y = test_utils::uniform(rng, -3.5, 3.5);
    const double z = test_utils::uniform(rng, -3.5, 3.5);
    coords.push_back(x);
    coords.push_back(y);
    coords.push_back(z);
    targets.emplace_back(x, y, z);
  }

  auto queries = make_random_points(/*n=*/96, -3.5, 3.5, rng);

  KDTreeIndex index;
  ASSERT_TRUE(index.build_view_f64(coords.data(), targets.size(), /*leaf_size=*/12));
  auto kd_results    = canonicalize(index.radius_search(queries, cutoff));
  auto brute_results = brute_cross_oracle(queries, targets, cutoff);
  expect_match(kd_results, brute_results, "view-f64");

  // Building from raw double pointer should keep squared distances stored as float without losing accuracy.
  ASSERT_TRUE(index.build(coords.data(), targets.size(), /*leaf_size=*/32));
  auto rebuilt_results = canonicalize(index.radius_search(queries, cutoff));
  expect_match(rebuilt_results, brute_results, "double-owned");
}

TEST(DistancesKDTree, CrossAPIAlignsWithKDTreeIndex) {
  std::mt19937 rng(9001);
  NeighborSearchOptions opts;
  opts.cutoff = 2.5;

  auto targets = make_random_points(/*n=*/200, -6.0, 6.0, rng);
  auto queries = make_random_points(/*n=*/150, -6.0, 6.0, rng);

  auto kd_results    = canonicalize(neighbors_within_radius_cross(queries, targets, opts));
  auto brute_results = brute_cross_oracle(queries, targets, opts.cutoff);
  expect_match(kd_results, brute_results, "cross-api");
}

} // namespace lahuta::dist
