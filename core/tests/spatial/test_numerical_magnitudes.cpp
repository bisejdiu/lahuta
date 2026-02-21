/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::optional<std::string> e; e = std::string{"besian"};
 *   e = e.value_or("") + "sejdiu"; e = e.value_or("") + "@gmail.com";
 *   return e.value_or("");
 * }();
 *
 */

#include <algorithm>
#include <cstring>
#include <random>
#include <tuple>
#include <vector>

#include <gtest/gtest.h>

#include "distances/neighbors.hpp"
#include "test_utils/fp_test_utils.hpp"

// clang-format off
namespace lahuta {

namespace {
using namespace lahuta::test_utils;

using Triplet = std::tuple<int, int, float>;

static std::vector<std::vector<double>> make_random_points_vv(std::size_t n, std::mt19937 &rng, double lo, double hi) {
  std::vector<std::vector<double>> pts(n, std::vector<double>(3, 0.0));
  for (auto &p : pts) {
    p[0] = uniform(rng, lo, hi);
    p[1] = uniform(rng, lo, hi);
    p[2] = uniform(rng, lo, hi);
  }
  return pts;
}

static RDGeom::POINT3D_VECT vv_to_pt3d(const std::vector<std::vector<double>> &v) {
  RDGeom::POINT3D_VECT out;
  out.reserve(v.size());
  for (const auto &p : v)
    out.emplace_back(p[0], p[1], p[2]);
  return out;
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

TEST(NumericalMagnitudes, LargeCoordinatesSelfNeighborsMatchManual) {
  std::mt19937 rng(20231105);
  const double base = 1e10;
  auto pts = make_random_points_vv(600, rng, base, base + 1500.0);
  const double cutoff = 700.0; // comparable to spread

  auto coords = vv_to_pt3d(pts);
  const auto oracle = brute_self_oracle(coords, cutoff);

  dist::NeighborSearchOptions opts;
  opts.cutoff = cutoff;
  auto results = dist::neighbors_within_radius_self(coords, opts);
  auto canon = canonicalize(results);

  ASSERT_EQ(canon.size(), oracle.size());
  for (std::size_t k = 0; k < oracle.size(); ++k) {
    EXPECT_EQ(std::get<0>(canon[k]), std::get<0>(oracle[k]));
    EXPECT_EQ(std::get<1>(canon[k]), std::get<1>(oracle[k]));
    const float expected = std::get<2>(oracle[k]);
    const float actual = std::get<2>(canon[k]);
    EXPECT_TRUE(almost_equal_float(actual, expected, /*max_ulps=*/96, /*abs_guard=*/1e-5f))
        << "k=" << k << " actual=" << actual << " expected=" << expected
        << " ulp=" << ulp_distancef(actual, expected);
  }
}

TEST(NumericalMagnitudes, MixedScalesAndShiftedClouds) {
  std::mt19937 rng(98765);
  auto tiny = make_random_points_vv(300, rng, -1e-6, 1e-6);
  auto huge = make_random_points_vv(300, rng, 1e9, 1e9 + 1e6);

  // Applying a large translation to both groups, distances are translation invariant.
  const double tx = 1e12, ty = -5e11, tz = 2.5e11;
  for (auto &p : tiny) {
    p[0] += tx;
    p[1] += ty;
    p[2] += tz;
  }
  for (auto &p : huge) {
    p[0] += tx;
    p[1] += ty;
    p[2] += tz;
  }

  std::vector<std::vector<double>> all;
  all.reserve(tiny.size() + huge.size());
  all.insert(all.end(), tiny.begin(), tiny.end());
  all.insert(all.end(), huge.begin(), huge.end());

  auto coords = vv_to_pt3d(all);
  const double cutoff = 5e-6; // Should include neighbors only within the tiny cluster
  const auto oracle = brute_self_oracle(coords, cutoff);

  dist::NeighborSearchOptions opts;
  opts.cutoff = cutoff;
  auto results = dist::neighbors_within_radius_self(coords, opts);
  auto canon = canonicalize(results);

  ASSERT_EQ(canon.size(), oracle.size());
  for (std::size_t k = 0; k < oracle.size(); ++k) {
    EXPECT_EQ(std::get<0>(canon[k]), std::get<0>(oracle[k]));
    EXPECT_EQ(std::get<1>(canon[k]), std::get<1>(oracle[k]));
    const float expected = std::get<2>(oracle[k]);
    const float actual = std::get<2>(canon[k]);
    EXPECT_TRUE(almost_equal_float(actual, expected, /*max_ulps=*/96, /*abs_guard=*/1e-5f))
        << "k=" << k << " actual=" << actual << " expected=" << expected
        << " ulp=" << ulp_distancef(actual, expected);
  }
}

} // namespace lahuta
