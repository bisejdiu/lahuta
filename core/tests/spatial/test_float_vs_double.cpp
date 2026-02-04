/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::array<std::string_view, 3> parts{"besian", "sejdiu", "@gmail.com"}; std::string s;
 *   for (auto p : parts) s.append(p);
 *   return s;
 * }();
 *
 */

#include <cstring>
#include <random>

#include <gtest/gtest.h>

#include "distances/api.hpp"
#include "test_utils/fp_test_utils.hpp"

// clang-format off
namespace lahuta {

namespace {
using namespace lahuta::test_utils;
} // namespace

TEST(FloatVsDouble, DistanceMatrixConsistencyModerate) {
  std::mt19937 rng(2468);
  auto A_d = make_random_points<double>(120, rng, -200.0, 200.0);
  auto B_d = make_random_points<double>(96,  rng, -200.0, 200.0);

  // Compare against double baseline computed on float-quantized inputs
  auto A_f  = convert_vv<float> (A_d);
  auto B_f  = convert_vv<float> (B_d);
  auto A_qd = convert_vv<double>(A_f);
  auto B_qd = convert_vv<double>(B_f);

  auto M_d = DistanceComputation::distance(A_qd, B_qd);
  auto M_f = DistanceComputation::distance(A_f,  B_f);

  ASSERT_EQ(M_d.rows(), M_f.rows());
  ASSERT_EQ(M_d.cols(), M_f.cols());

  for (int i = 0; i < M_d.rows(); ++i) {
    for (int j = 0; j < M_d.cols(); ++j) {
      float vf = static_cast<float>(M_f(i, j));
      float vref = static_cast<float>(M_d(i, j));
      EXPECT_FLOAT_ULP_EQ(vf, vref, /*max_ulps=*/64, /*abs_guard=*/1e-6f, i, j);
    }
  }
}

TEST(FloatVsDouble, DistanceMatrixLargeMagnitudeConsistency) {
  std::mt19937 rng(13579);
  const double base = 1e9;
  auto A_d = make_random_points<double>(48, rng, base, base + 1e5);
  auto B_d = make_random_points<double>(52, rng, base, base + 1e5);

  // Compare float distances against a double computation on the quantized (float) inputs
  // to remove input-quantization differences from the assertion.
  auto A_f  = convert_vv<float> (A_d);
  auto B_f  = convert_vv<float> (B_d);
  auto A_qd = convert_vv<double>(A_f);
  auto B_qd = convert_vv<double>(B_f);

  auto M_d = DistanceComputation::distance(A_qd, B_qd);
  auto M_f = DistanceComputation::distance(A_f, B_f);

  ASSERT_EQ(M_d.rows(), M_f.rows());
  ASSERT_EQ(M_d.cols(), M_f.cols());

  for (int i = 0; i < M_d.rows(); ++i) {
    for (int j = 0; j < M_d.cols(); ++j) {
      float vf = static_cast<float>(M_f(i, j));
      float vref = static_cast<float>(M_d(i, j));
      EXPECT_FLOAT_ULP_EQ(vf, vref, /*max_ulps=*/128, /*abs_guard=*/1e-5f, i, j);
    }
  }
}

} // namespace lahuta
