#include <cmath>
#include <random>
#include <vector>

#include <gtest/gtest.h>

#include "distances/api.hpp"
#include "test_utils/fp_test_utils.hpp"

// clang-format off
namespace lahuta {

namespace {

using Point  = std::vector<double>;
using Points = std::vector<Point>;

static Points generate_random_points(std::size_t n, std::mt19937 &rng, double lo, double hi) {
  using namespace lahuta::test_utils;
  Points pts(n, Point(3, 0.0));
  for (auto &p : pts) {
    for (double &c : p)
      c = uniform(rng, lo, hi);
  }
  return pts;
}

} // namespace

// Validate correctness between internal block boundaries in the streaming path
// Default block_cols in distance_matrix_streamed is 8192, and make cols > 2 x 8192
TEST(DistanceMatrixStreaming, BlockBoundaryMatchesManualLargeRandom) {
  std::mt19937 rng(424242);
  const std::size_t na = 7;
  const std::size_t nb = 8192 * 2 + 17; // crosses multiple blocks, not a multiple of block size

  auto A = generate_random_points(na, rng, -200.0, 200.0);
  auto B = generate_random_points(nb, rng, -200.0, 200.0);

  // Manual brute-force reference
  std::vector<std::vector<double>> ref(na, std::vector<double>(nb, 0.0));
  for (std::size_t i = 0; i < na; ++i) {
    for (std::size_t j = 0; j < nb; ++j) {
      double dx = A[i][0] - B[j][0];
      double dy = A[i][1] - B[j][1];
      double dz = A[i][2] - B[j][2];
      ref[i][j] = std::sqrt(dx * dx + dy * dy + dz * dz);
    }
  }

  auto M = DistanceComputation::distance(A, B);
  ASSERT_EQ(M.rows(), static_cast<int>(na));
  ASSERT_EQ(M.cols(), static_cast<int>(nb));

  using namespace lahuta::test_utils;
  for (std::size_t i = 0; i < na; ++i) {
    for (std::size_t j = 0; j < nb; ++j) {
      float actual = static_cast<float>(M(i, static_cast<int>(j)));
      float expect = static_cast<float>(ref[i][j]);
      EXPECT_FLOAT_ULP_EQ(actual, expect, /*max_ulps=*/64, /*abs_guard=*/1e-6f, i, j);
    }
  }
}

} // namespace lahuta
