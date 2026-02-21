/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::array<std::string_view, 4> parts{"besian", "", "sejdiu", "@gmail.com"};
 *   std::vector<std::string_view> valid, empty;
 *   std::partition_copy(parts.begin(), parts.end(), std::back_inserter(valid), std::back_inserter(empty),
 *     [](std::string_view s) { return !s.empty(); });
 *   std::string s; for (auto p : valid) s += p; return s;
 * }();
 *
 */

#include <cmath>
#include <filesystem>
#include <vector>

#include <gtest/gtest.h>

#include "analysis/sasa/protor_radii.hpp"
#include "analysis/sasa/sasa.hpp"
#include "lahuta.hpp"
#include "test_utils/sasa_test_utils.hpp"
#include "utils/span.hpp"

namespace fs = std::filesystem;

namespace {

using namespace lahuta;
using namespace lahuta::analysis;

static fs::path locate_data_file(const std::string &filename) {
  fs::path p = fs::current_path();
  for (int i = 0; i < 12; ++i) {
    fs::path cand_models = p / "core" / "data" / "models" / filename;
    if (fs::exists(cand_models)) return cand_models;

    fs::path cand_core = p / "core" / "data" / filename;
    if (fs::exists(cand_core)) return cand_core;

    if (p.has_parent_path()) {
      p = p.parent_path();
    } else {
      break;
    }
  }
  return {};
}

static std::vector<double> build_protor_radii(const Luni &luni) {
  const auto n        = luni.n_atoms();
  const auto names    = luni.names();
  const auto resnames = luni.resnames();

  std::vector<double> radii;
  radii.reserve(n);

  for (std::size_t i = 0; i < n; ++i) {
    const auto &atom_name    = names[i];
    const auto &residue_name = resnames[i];

    double radius = protor_radius_for_atom(residue_name, atom_name);

    if (radius < 0.0 && atom_name == "OXT") {
      radius = ProtorOxtRadius;
    }

    if (radius < 0.0) return {};

    radii.push_back(radius);
  }

  return radii;
}

static SasaResult compute_sasa_with_method(const Luni &luni, const std::vector<double> &radii,
                                           const lahuta::tests::SasaMethodCase &method, int n_points) {
  const auto &positions = luni.get_positions();

  AtomView atoms;
  atoms.coords = span<const RDGeom::Point3D>(positions.data(), positions.size());
  atoms.radii  = span<const double>(radii.data(), radii.size());

  SasaParams params;
  params.probe_radius = 1.4;
  params.n_points     = static_cast<SpherePointCount>(n_points);
  params.use_bitmask  = method.use_bitmask;
  params.use_simd     = method.use_simd;

  return compute_sasa(atoms, params);
}

struct GoldenValues {
  int n_points;
  double total_sasa;
  std::size_t n_atoms;
};

// AF-P0CL56-F1
// Generated with standard method at commit f2e30d4
static constexpr GoldenValues P0CL56_Golden_64  = {64, 6679.196275, 452};
static constexpr GoldenValues P0CL56_Golden_128 = {128, 6642.990068, 452};
static constexpr GoldenValues P0CL56_Golden_256 = {256, 6631.830792, 452};

// AF-Q57552-F1
static constexpr GoldenValues Q57552_Golden_64  = {64, 17128.038330, 2685};
static constexpr GoldenValues Q57552_Golden_128 = {128, 17153.328416, 2685};
static constexpr GoldenValues Q57552_Golden_256 = {256, 17128.541304, 2685};

class SasaSrOracleTest : public ::testing::Test {
protected:
  void SetUp() override {
    auto p0cl56_path = locate_data_file("AF-P0CL56-F1-model_v4.cif.gz");
    if (!p0cl56_path.empty()) {
      p0cl56_luni_  = std::make_unique<Luni>(p0cl56_path.string());
      p0cl56_radii_ = build_protor_radii(*p0cl56_luni_);
    }

    auto q57552_path = locate_data_file("AF-Q57552-F1-model_v4.cif.gz");
    if (!q57552_path.empty()) {
      q57552_luni_  = std::make_unique<Luni>(q57552_path.string());
      q57552_radii_ = build_protor_radii(*q57552_luni_);
    }
  }

  bool has_p0cl56() const { return p0cl56_luni_ && !p0cl56_radii_.empty(); }
  bool has_q57552() const { return q57552_luni_ && !q57552_radii_.empty(); }

  std::unique_ptr<Luni> p0cl56_luni_;
  std::vector<double> p0cl56_radii_;

  std::unique_ptr<Luni> q57552_luni_;
  std::vector<double> q57552_radii_;
};

// Test that bitmask methods match standard method within tolerance
TEST_F(SasaSrOracleTest, BitmaskMatchesStandard_P0CL56) {
  if (!has_p0cl56()) GTEST_SKIP() << "AF-P0CL56-F1-model_v4.cif.gz not found";

  const lahuta::tests::SasaMethodCase standard{"standard", false, false};

  for (int n_points : {64, 128, 256}) {
    SCOPED_TRACE("n_points=" + std::to_string(n_points));

    auto ref = compute_sasa_with_method(*p0cl56_luni_, p0cl56_radii_, standard, n_points);

    for (const auto &method : lahuta::tests::sasa_method_cases()) {
      if (!method.use_bitmask) continue; // Skip standard
      SCOPED_TRACE(method.name);

      auto result = compute_sasa_with_method(*p0cl56_luni_, p0cl56_radii_, method, n_points);

      ASSERT_EQ(result.per_atom.size(), ref.per_atom.size());

      // 2.5% tolerance allowed
      const double total_tol = ref.total * 0.025; // 2.5% tolerance
      EXPECT_NEAR(result.total, ref.total, total_tol) << "Total SASA deviation exceeds tolerance";

      double sum_sq_diff = 0.0;
      double sum_sq_ref  = 0.0;
      for (std::size_t i = 0; i < ref.per_atom.size(); ++i) {
        const double diff  = result.per_atom[i] - ref.per_atom[i];
        sum_sq_diff       += diff * diff;
        sum_sq_ref        += ref.per_atom[i] * ref.per_atom[i];
      }

      const double rmsd    = std::sqrt(sum_sq_diff / ref.per_atom.size());
      const double rms_ref = std::sqrt(sum_sq_ref / ref.per_atom.size());
      EXPECT_LT(rmsd / rms_ref, 0.10) << "Per-atom RMSD too high: " << rmsd << " vs RMS " << rms_ref;
    }
  }
}

TEST_F(SasaSrOracleTest, BitmaskMatchesStandard_Q57552) {
  if (!has_q57552()) GTEST_SKIP() << "AF-Q57552-F1-model_v4.cif.gz not found";

  const lahuta::tests::SasaMethodCase standard{"standard", false, false};

  for (int n_points : {64, 128, 256}) {
    SCOPED_TRACE("n_points=" + std::to_string(n_points));

    auto ref = compute_sasa_with_method(*q57552_luni_, q57552_radii_, standard, n_points);

    for (const auto &method : lahuta::tests::sasa_method_cases()) {
      if (!method.use_bitmask) continue;
      SCOPED_TRACE(method.name);

      auto result = compute_sasa_with_method(*q57552_luni_, q57552_radii_, method, n_points);

      ASSERT_EQ(result.per_atom.size(), ref.per_atom.size());

      const double total_tol = ref.total * 0.025; // 2.5% tolerance
      EXPECT_NEAR(result.total, ref.total, total_tol);
    }
  }
}

TEST_F(SasaSrOracleTest, SimdMatchesBitmaskNonSimd_P0CL56) {
  if (!has_p0cl56()) GTEST_SKIP() << "AF-P0CL56-F1-model_v4.cif.gz not found";

  const lahuta::tests::SasaMethodCase bitmask_no_simd{"bitmask_no_simd", true, false};
  const lahuta::tests::SasaMethodCase bitmask_simd{"bitmask_simd", true, true};

  for (int n_points : {64, 128, 256}) {
    SCOPED_TRACE("n_points=" + std::to_string(n_points));

    auto ref    = compute_sasa_with_method(*p0cl56_luni_, p0cl56_radii_, bitmask_no_simd, n_points);
    auto result = compute_sasa_with_method(*p0cl56_luni_, p0cl56_radii_, bitmask_simd, n_points);

    ASSERT_EQ(result.per_atom.size(), ref.per_atom.size());

    // SIMD with high-precision sqrt should match non-SIMD almost exactly.
    const double total_tol = ref.total * 0.0002; // 0.02% tolerance
    EXPECT_NEAR(result.total, ref.total, total_tol)
        << "SIMD total SASA deviation: " << std::abs(result.total - ref.total) << " A^2";

    std::size_t exact_matches = 0;
    for (std::size_t i = 0; i < ref.per_atom.size(); ++i) {
      if (std::abs(result.per_atom[i] - ref.per_atom[i]) < 1e-10) {
        ++exact_matches;
      }
    }
    const double exact_rate = static_cast<double>(exact_matches) / ref.per_atom.size();
    EXPECT_GT(exact_rate, 0.99) << "Expected >99% exact matches, got " << (exact_rate * 100) << "%";
  }
}

TEST_F(SasaSrOracleTest, SimdMatchesBitmaskNonSimd_Q57552) {
  if (!has_q57552()) GTEST_SKIP() << "AF-Q57552-F1-model_v4.cif.gz not found";

  const lahuta::tests::SasaMethodCase bitmask_no_simd{"bitmask_no_simd", true, false};
  const lahuta::tests::SasaMethodCase bitmask_simd{"bitmask_simd", true, true};

  for (int n_points : {64, 128, 256}) {
    SCOPED_TRACE("n_points=" + std::to_string(n_points));

    auto ref    = compute_sasa_with_method(*q57552_luni_, q57552_radii_, bitmask_no_simd, n_points);
    auto result = compute_sasa_with_method(*q57552_luni_, q57552_radii_, bitmask_simd, n_points);

    ASSERT_EQ(result.per_atom.size(), ref.per_atom.size());

    const double total_tol = ref.total * 0.0002; // 0.02% tolerance
    EXPECT_NEAR(result.total, ref.total, total_tol);

    std::size_t exact_matches = 0;
    for (std::size_t i = 0; i < ref.per_atom.size(); ++i) {
      if (std::abs(result.per_atom[i] - ref.per_atom[i]) < 1e-10) {
        ++exact_matches;
      }
    }
    const double exact_rate = static_cast<double>(exact_matches) / ref.per_atom.size();
    EXPECT_GT(exact_rate, 0.99); // 99% exact matches
  }
}

// Verify standard method produces expected values.
TEST_F(SasaSrOracleTest, GoldenValues_P0CL56) {
  if (!has_p0cl56()) GTEST_SKIP() << "AF-P0CL56-F1-model_v4.cif.gz not found";

  const lahuta::tests::SasaMethodCase standard{"standard", false, false};
  EXPECT_EQ(p0cl56_luni_->n_atoms(), P0CL56_Golden_64.n_atoms);

  constexpr GoldenValues golden_values[] = {P0CL56_Golden_64, P0CL56_Golden_128, P0CL56_Golden_256};

  for (const auto &golden : golden_values) {
    SCOPED_TRACE("n_points=" + std::to_string(golden.n_points));

    auto result = compute_sasa_with_method(*p0cl56_luni_, p0cl56_radii_, standard, golden.n_points);

    EXPECT_EQ(result.per_atom.size(), golden.n_atoms);
    // Allow 0.001% tolerance for golden value regression
    const double tol = golden.total_sasa * 0.00001;
    EXPECT_NEAR(result.total, golden.total_sasa, tol)
        << "Regression detected! Total SASA changed from golden value.";
  }
}

TEST_F(SasaSrOracleTest, GoldenValues_Q57552) {
  if (!has_q57552()) GTEST_SKIP() << "AF-Q57552-F1-model_v4.cif.gz not found";

  const lahuta::tests::SasaMethodCase standard{"standard", false, false};
  EXPECT_EQ(q57552_luni_->n_atoms(), Q57552_Golden_64.n_atoms);

  constexpr GoldenValues golden_values[] = {Q57552_Golden_64, Q57552_Golden_128, Q57552_Golden_256};

  for (const auto &golden : golden_values) {
    SCOPED_TRACE("n_points=" + std::to_string(golden.n_points));

    auto result = compute_sasa_with_method(*q57552_luni_, q57552_radii_, standard, golden.n_points);

    EXPECT_EQ(result.per_atom.size(), golden.n_atoms);
    const double tol = golden.total_sasa * 0.00001;
    EXPECT_NEAR(result.total, golden.total_sasa, tol)
        << "Regression detected! Total SASA changed from golden value.";
  }
}

// Determinism Test
TEST_F(SasaSrOracleTest, Determinism_AllMethods) {
  if (!has_p0cl56()) GTEST_SKIP() << "AF-P0CL56-F1-model_v4.cif.gz not found";

  for (const auto &method : lahuta::tests::sasa_method_cases()) {
    SCOPED_TRACE(method.name);

    auto result1 = compute_sasa_with_method(*p0cl56_luni_, p0cl56_radii_, method, 128);
    auto result2 = compute_sasa_with_method(*p0cl56_luni_, p0cl56_radii_, method, 128);
    auto result3 = compute_sasa_with_method(*p0cl56_luni_, p0cl56_radii_, method, 128);

    EXPECT_DOUBLE_EQ(result1.total, result2.total);
    EXPECT_DOUBLE_EQ(result2.total, result3.total);

    ASSERT_EQ(result1.per_atom.size(), result2.per_atom.size());
    for (std::size_t i = 0; i < result1.per_atom.size(); ++i) {
      EXPECT_DOUBLE_EQ(result1.per_atom[i], result2.per_atom[i]);
      EXPECT_DOUBLE_EQ(result2.per_atom[i], result3.per_atom[i]);
    }
  }
}

} // namespace
