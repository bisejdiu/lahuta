#include <gtest/gtest.h>

#include <algorithm>
#include <random>
#include <tuple>
#include <vector>

#include "distances/neighbors.hpp"
#include "test_utils/fp_test_utils.hpp"

// clang-format off
namespace lahuta::dist {

namespace {

using Triplet = std::tuple<int, int, float>; // (query_idx, target_idx, d2)

static RDGeom::POINT3D_VECT make_degenerate_targets(int n, double spacing) {
  RDGeom::POINT3D_VECT t;
  t.reserve(n);
  for (int i = 0; i < n; ++i) {
    t.emplace_back(i * spacing, 0.0, 0.0); // all along x-axis
  }
  return t;
}

static RDGeom::POINT3D_VECT make_thin_queries(int n, double x_min, double x_max, double yz_extent, std::mt19937 &rng) {
  using namespace lahuta::test_utils;
  RDGeom::POINT3D_VECT q;
  q.reserve(n);
  for (int i = 0; i < n; ++i) {
    q.emplace_back(
        uniform(rng, x_min, x_max),
        uniform(rng, -yz_extent, yz_extent),
        uniform(rng, -yz_extent, yz_extent));
  }
  return q;
}

static std::vector<Triplet> brute_cross_oracle(const RDGeom::POINT3D_VECT &queries, const RDGeom::POINT3D_VECT &targets, double cutoff) {
  const float cutoff_sq = static_cast<float>(cutoff * cutoff);
  std::vector<Triplet> out;
  const int nq = static_cast<int>(queries.size());
  const int nt = static_cast<int>(targets.size());
  out.reserve(std::min(nq * nt, 200000));
  for (int i = 0; i < nq; ++i) {
    for (int j = 0; j < nt; ++j) {
      double dx = queries[i].x - targets[j].x;
      double dy = queries[i].y - targets[j].y;
      double dz = queries[i].z - targets[j].z;
      float d2 = static_cast<float>(dx * dx + dy * dy + dz * dz);
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
    out.emplace_back(pairs[k].first, pairs[k].second, dists[k]);
  }
  std::sort(out.begin(), out.end());
  return out;
}

} // namespace

// Force cross-fallback by making targets degenerate (1D line), but queries thin 3D cloud.
TEST(CrossNeighborsFallbackStrict, TargetsDegenerateQueriesThinCloud) {
  std::mt19937 rng(424242);
  const int nt = 256;
  const int nq = 220;
  const double spacing = 0.25;    // target spacing along x
  const double x_min   = 0.0, x_max = spacing * (nt - 1);
  const double yz_extent = 0.10;  // queries hover near x-axis
  const double cutoff = 0.55;     // captures nearby targets even with small yz offsets

  auto targets = make_degenerate_targets(nt, spacing);
  auto queries = make_thin_queries(nq, x_min, x_max, yz_extent, rng);

  NeighborSearchOptions opts;
  opts.cutoff = cutoff;

  // Build will fail for targets due to degenerate dimensions, forcing brute-force cross fallback.
  auto results = neighbors_within_radius_cross(queries, targets, opts);
  auto canon   = canonicalize(results);
  auto oracle  = brute_cross_oracle(queries, targets, cutoff);

  ASSERT_EQ(canon.size(), oracle.size());
  for (std::size_t k = 0; k < oracle.size(); ++k) {
    EXPECT_EQ(std::get<0>(canon[k]), std::get<0>(oracle[k]));
    EXPECT_EQ(std::get<1>(canon[k]), std::get<1>(oracle[k]));
    EXPECT_FLOAT_ULP_EQ(
      std::get<2>(canon[k]),
      std::get<2>(oracle[k]),
      /*max_ulps=*/256,
      /*abs_guard=*/1e-6f,
      std::get<0>(oracle[k]),
      std::get<1>(oracle[k])
    );
  }
}

} // namespace lahuta::dist
