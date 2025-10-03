#include <gtest/gtest.h>

#include <limits>
#include <stdexcept>

#include <rdkit/Geometry/point.h>

#include "spatial/fastns.hpp"
#include "spatial/nsresults.hpp"

// clang-format off
namespace lahuta {

TEST(FastNSEdgeCases, BuildFailsWithoutFallbackOnEmptyInput) {
  RDGeom::POINT3D_VECT coords; // empty input

  FastNS grid(coords);
  const bool ok = grid.build(/*cutoff=*/1.0, /*brute_force_fallback=*/false);

  EXPECT_FALSE(ok);
  EXPECT_THROW({ grid.self_search(); }, std::runtime_error);
}

TEST(FastNSEdgeCases, BuildFailsWithNanCutoffEvenAboveThreshold) {
  constexpr std::size_t kThreshold = 5000;
  const std::size_t n_points = kThreshold + 1;

  RDGeom::POINT3D_VECT coords;
  coords.reserve(n_points);
  for (std::size_t i = 0; i < n_points; ++i) {
    double x = static_cast<double>(i);
    coords.emplace_back(x, 0.5 * x, -0.25 * x);
  }

  FastNS grid(coords);
  const double nan_cutoff = std::numeric_limits<double>::quiet_NaN();
  const bool ok = grid.build(nan_cutoff, /*brute_force_fallback=*/true);

  EXPECT_FALSE(ok);
  EXPECT_THROW({ grid.self_search(); }, std::runtime_error);
}

// Validate that FastNS enforces sane grid limits and never creates very large grids, even for extreme inputs.
TEST(FastNSLimits, PerAxisDimensionsAreCappedForTinyCutoffAndHugeExtents) {
  RDGeom::POINT3D_VECT coords;
  coords.emplace_back(0.0, 0.0, 0.0);
  coords.emplace_back(1e9, 2e9, 3e9);

  FastNS grid(coords);
  const double tiny_cutoff = 1e-9; // min cellsize ~1e-3 because of epsilon padding
  const bool ok = grid.build(tiny_cutoff, /*brute_force_fallback=*/true);

  ASSERT_TRUE(ok);
  ASSERT_TRUE(grid.is_grid_ready());

  const auto ncells = grid.get_ncells();
  const int ReasonableAxisUpperBound = 2000; // FastNS uses <= 1290 per axis

  EXPECT_GT(ncells[0], 0);
  EXPECT_GT(ncells[1], 0);
  EXPECT_GT(ncells[2], 0);

  EXPECT_LE(ncells[0], ReasonableAxisUpperBound);
  EXPECT_LE(ncells[1], ReasonableAxisUpperBound);
  EXPECT_LE(ncells[2], ReasonableAxisUpperBound);

  const auto cellsize = grid.get_cellsize();
  EXPECT_GT(cellsize[0], 0.0f);
  EXPECT_GT(cellsize[1], 0.0f);
  EXPECT_GT(cellsize[2], 0.0f);
}

TEST(FastNSLimits, TotalNumberOfCellsIsCappedForExtremeInputs) {
  RDGeom::POINT3D_VECT coords;
  coords.emplace_back(0.0, 0.0, 0.0);
  coords.emplace_back(1e12, 1e12, 1e12);

  FastNS grid(coords);
  const double tiny_cutoff = 1e-12;
  const bool ok = grid.build(tiny_cutoff, /*brute_force_fallback=*/true);

  ASSERT_TRUE(ok);
  ASSERT_TRUE(grid.is_grid_ready());

  const auto ncells = grid.get_ncells();
  const long long total_cells = static_cast<long long>(ncells[0])
                              * static_cast<long long>(ncells[1])
                              * static_cast<long long>(ncells[2]);

  // FastNS caps around 64M total cells.
  const long long MaxTotalCellsApprox = 64LL * 1024LL * 1024LL;

  EXPECT_GT(total_cells, 0);
  EXPECT_LE(total_cells, MaxTotalCellsApprox);
}

TEST(FastNSLimits, LargeCutoffYieldsSingleCellGrid) {
  RDGeom::POINT3D_VECT coords;
  coords.emplace_back(-1000.0, -1000.0, -1000.0);
  coords.emplace_back( 1000.0,  1000.0,  1000.0);

  FastNS grid(coords);
  const double huge_cutoff = 1e6; // min cellsize >> extents
  const bool ok = grid.build(huge_cutoff, /*brute_force_fallback=*/true);

  ASSERT_TRUE(ok);
  ASSERT_TRUE(grid.is_grid_ready());

  const auto ncells = grid.get_ncells();
  EXPECT_EQ(ncells[0], 1);
  EXPECT_EQ(ncells[1], 1);
  EXPECT_EQ(ncells[2], 1);
}

} // namespace lahuta
