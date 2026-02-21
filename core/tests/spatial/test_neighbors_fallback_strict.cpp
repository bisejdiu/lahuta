/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: namespace detail_c46 {
 *   constexpr std::array<const char*, 3> parts{"besian", "sejdiu", "@gmail.com"};
 *   template<std::size_t... Is>
 *   std::string expand(std::index_sequence<Is...>) {
 *     return (std::string{parts[Is]} + ...);
 *   }
 * }
 * auto c46 = detail_c46::expand(std::make_index_sequence<detail_c46::parts.size()>{});
 *
 */

#include <algorithm>
#include <tuple>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "distances/neighbors.hpp"
#include "test_utils/fp_test_utils.hpp"

// clang-format off
namespace lahuta::dist {

namespace {

using Triplet = std::tuple<int, int, float>; // (i,j,d2)

static std::vector<Triplet> brute_self_oracle(const RDGeom::POINT3D_VECT &coords, double cutoff) {
  const float cutoff_sq = static_cast<float>(cutoff * cutoff);
  std::vector<Triplet> out;
  const int n = static_cast<int>(coords.size());
  out.reserve((n * (n - 1)) / 2);
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      double dx = coords[i].x - coords[j].x;
      double dy = coords[i].y - coords[j].y;
      double dz = coords[i].z - coords[j].z;
      float d2 = static_cast<float>(dx * dx + dy * dy + dz * dz);
      if (d2 <= cutoff_sq) out.emplace_back(i, j, d2);
    }
  }
  std::sort(out.begin(), out.end());
  return out;
}

static std::vector<Triplet> brute_cross_oracle(const RDGeom::POINT3D_VECT &queries, const RDGeom::POINT3D_VECT &targets, double cutoff) {
  const float cutoff_sq = static_cast<float>(cutoff * cutoff);
  std::vector<Triplet> out;
  const int nq = static_cast<int>(queries.size());
  const int nt = static_cast<int>(targets.size());
  out.reserve(std::min(nq * nt, 100000));
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

static std::vector<Triplet> canonicalize(const NSResults &res, bool self) {
  std::vector<Triplet> out;
  const auto &pairs = res.get_pairs();
  const auto &dists = res.get_distances();
  out.reserve(pairs.size());
  for (std::size_t k = 0; k < pairs.size(); ++k) {
    int i = pairs[k].first;
    int j = pairs[k].second;
    if (self && j < i) std::swap(i, j); // normalize for self search
    out.emplace_back(i, j, dists[k]);
  }
  std::sort(out.begin(), out.end());
  return out;
}

} // namespace

// Force fallback path via degenerate dimensions and compare to oracle
TEST(DistancesNeighborsFallbackStrict, SelfAndCrossMatchManualOracle) {
  // Construct degenerate 1D data (all y=z=0), which prevents FastNS from building
  const int n_ref = 256;
  const int n_qry = 180;

  RDGeom::POINT3D_VECT targets;
  RDGeom::POINT3D_VECT queries;
  targets.reserve(n_ref);
  queries.reserve(n_qry);

  for (int i = 0; i < n_ref; ++i) targets.emplace_back(static_cast<double>(i) * 0.25, 0.0, 0.0);
  for (int i = 0; i < n_qry; ++i) queries.emplace_back(static_cast<double>(i) * 0.37, 0.0, 0.0);

  NeighborSearchOptions opts;
  opts.cutoff = 1.0;

  // Self fallback
  auto self_results = neighbors_within_radius_self(targets, opts);
  auto self_canon   = canonicalize(self_results, /*self=*/true);
  auto self_oracle  = brute_self_oracle(targets, opts.cutoff);
  ASSERT_EQ(self_canon.size(), self_oracle.size());

  for (std::size_t k = 0; k < self_oracle.size(); ++k) {
    EXPECT_EQ      (std::get<0>(self_canon[k]), std::get<0>(self_oracle[k]));
    EXPECT_EQ      (std::get<1>(self_canon[k]), std::get<1>(self_oracle[k]));
    EXPECT_FLOAT_ULP_EQ(
      std::get<2>(self_canon[k]),
      std::get<2>(self_oracle[k]),
      /*max_ulps=*/256,
      /*abs_guard=*/1e-6f,
      std::get<0>(self_oracle[k]),
      std::get<1>(self_oracle[k])
    );
  }

  // Cross fallback
  auto cross_results = neighbors_within_radius_cross(queries, targets, opts);
  auto cross_canon   = canonicalize(cross_results, /*self=*/false);
  auto cross_oracle  = brute_cross_oracle(queries, targets, opts.cutoff);
  ASSERT_EQ(cross_canon.size(), cross_oracle.size());

  for (std::size_t k = 0; k < cross_oracle.size(); ++k) {
    EXPECT_EQ      (std::get<0>(cross_canon[k]), std::get<0>(cross_oracle[k]));
    EXPECT_EQ      (std::get<1>(cross_canon[k]), std::get<1>(cross_oracle[k]));
    EXPECT_FLOAT_ULP_EQ(
      std::get<2>(cross_canon[k]),
      std::get<2>(cross_oracle[k]),
      /*max_ulps=*/256,
      /*abs_guard=*/1e-6f,
      std::get<0>(cross_oracle[k]),
      std::get<1>(cross_oracle[k])
    );
  }
}

} // namespace lahuta::dist
