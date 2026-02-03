/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto f = [](auto arg) {
 *     using T = std::decay_t<decltype(arg)>;
 *     if constexpr (std::disjunction_v<std::is_same<T, const char*>, std::is_same<T, std::string_view>>) return std::string(arg);
 *   };
 *   return f("besian") + f("sejdiu") + f("@gmail.com");
 * }();
 *
 */

#include <cmath>
#include <vector>

#include <Geometry/point.h>
#include <gtest/gtest.h>

#include "analysis/sasa/sasa.hpp"
#include "utils/math_constants.hpp"
#include "utils/span.hpp"
#include "test_utils/sasa_test_utils.hpp"

namespace {

using lahuta::Pi;

double sphere_area(double r) { return 4.0 * Pi * r * r; }

double expected_overlap_area_equal_spheres(double r, double center_dist) {
  if (center_dist >= 2.0 * r) return sphere_area(r);
  if (center_dist <= 0.0) return 0.0;
  const double cap_height = r - 0.5 * center_dist;
  const double occluded   = 2.0 * Pi * r * cap_height;
  return sphere_area(r) - occluded;
}

std::vector<double> quantize_radii(const std::vector<double> &radii) {
  std::vector<double> out;
  out.reserve(radii.size());
  for (double v : radii) {
    out.push_back(static_cast<double>(static_cast<float>(v)));
  }
  return out;
}

std::vector<RDGeom::Point3D> quantize_points(const std::vector<RDGeom::Point3D> &coords) {
  std::vector<RDGeom::Point3D> out;
  out.reserve(coords.size());
  for (const auto &pt : coords) {
    const float xf = static_cast<float>(pt.x);
    const float yf = static_cast<float>(pt.y);
    const float zf = static_cast<float>(pt.z);
    out.emplace_back(static_cast<double>(xf), static_cast<double>(yf), static_cast<double>(zf));
  }
  return out;
}

lahuta::analysis::AtomView make_atoms(const std::vector<RDGeom::Point3D> &coords,
                                      const std::vector<double> &radii) {
  return {lahuta::span<const RDGeom::Point3D>(coords.data(), coords.size()),
          lahuta::span<const double>(radii.data(), radii.size())};
}

} // namespace

TEST(SasaSrPrecision, FloatQuantizedInputsMatchWithinTolerance) {
  std::vector<RDGeom::Point3D> coords = {
      RDGeom::Point3D(0.1, 0.2, 0.3),
      RDGeom::Point3D(3.25, 0.1, -0.2),
      RDGeom::Point3D(0.2, 4.15, 0.05),
  };
  std::vector<double> radii = {2.0, 2.0, 2.0};

  const auto coords_q = quantize_points(coords);
  const auto radii_q  = quantize_radii(radii);

  lahuta::analysis::SasaParams base;
  base.probe_radius = 0.5;
  base.n_points     = 256;

  const auto atoms_d  = make_atoms(coords, radii);
  const auto atoms_qd = make_atoms(coords_q, radii_q);

  for (const auto &method : lahuta::tests::sasa_method_cases()) {
    SCOPED_TRACE(method.name);
    const auto params   = lahuta::tests::apply_sasa_method(base, method);
    const auto result_d = lahuta::analysis::compute_sasa(atoms_d, params);
    const auto result_qd = lahuta::analysis::compute_sasa(atoms_qd, params);

    ASSERT_EQ(result_d.per_atom.size(), result_qd.per_atom.size());

    constexpr double rel_tol = 0.01;
    double total_tol         = 0.0;
    for (std::size_t i = 0; i < result_d.per_atom.size(); ++i) {
      const double full_area = sphere_area(radii[i] + params.probe_radius);
      const double tol       = full_area * rel_tol;
      total_tol             += tol;
      EXPECT_NEAR(result_qd.per_atom[i], result_d.per_atom[i], tol);
    }

    EXPECT_NEAR(result_qd.total, result_d.total, total_tol);
  }
}

TEST(SasaSrDeterminism, SameOutputAcrossRuns) {
  std::vector<RDGeom::Point3D> coords = {
      RDGeom::Point3D(0.0, 0.0, 0.0),
      RDGeom::Point3D(3.0, 0.0, 0.0),
      RDGeom::Point3D(0.0, 4.5, 0.0),
  };
  std::vector<double> radii = {2.0, 2.0, 1.5};

  const auto atoms = make_atoms(coords, radii);

  lahuta::analysis::SasaParams base;
  base.probe_radius = 0.5;
  base.n_points     = 256;

  for (const auto &method : lahuta::tests::sasa_method_cases()) {
    SCOPED_TRACE(method.name);
    const auto params = lahuta::tests::apply_sasa_method(base, method);
    const auto ref    = lahuta::analysis::compute_sasa(atoms, params);

    for (int iter = 0; iter < 3; ++iter) {
      const auto res = lahuta::analysis::compute_sasa(atoms, params);
      ASSERT_EQ(res.per_atom.size(), ref.per_atom.size());
      for (std::size_t i = 0; i < res.per_atom.size(); ++i) {
        EXPECT_NEAR(res.per_atom[i], ref.per_atom[i], 1e-12);
      }
      EXPECT_NEAR(res.total, ref.total, 1e-12);
    }
  }
}

TEST(SasaSrRegression, FixedStructureMatchesGoldenOutputs) {
  std::vector<RDGeom::Point3D> coords = {
      RDGeom::Point3D(0.0, 0.0, 0.0),
      RDGeom::Point3D(3.0, 0.0, 0.0),
      RDGeom::Point3D(10.0, 0.0, 0.0),
      RDGeom::Point3D(0.0, 8.0, 0.0),
  };
  std::vector<double> radii = {2.0, 2.0, 1.5, 1.0};

  lahuta::analysis::SasaParams base;
  base.probe_radius = 0.5;
  base.n_points     = 256;

  const auto atoms = make_atoms(coords, radii);

  for (const auto &method : lahuta::tests::sasa_method_cases()) {
    SCOPED_TRACE(method.name);
    const auto params  = lahuta::tests::apply_sasa_method(base, method);
    const auto result  = lahuta::analysis::compute_sasa(atoms, params);

    const double r_overlap        = 2.0 + params.probe_radius;
    const double expected_overlap = expected_overlap_area_equal_spheres(r_overlap, 3.0);
    const double expected_far_1   = sphere_area(1.5 + params.probe_radius);
    const double expected_far_2   = sphere_area(1.0 + params.probe_radius);

    const std::vector<double> expected = {
        expected_overlap,
        expected_overlap,
        expected_far_1,
        expected_far_2,
    };

    ASSERT_EQ(result.per_atom.size(), expected.size());

    constexpr double rel_tol = 0.05;
    double expected_total    = 0.0;
    for (std::size_t i = 0; i < expected.size(); ++i) {
      expected_total      += expected[i];
      const double denom   = std::max(1e-12, std::abs(expected[i]));
      const double rel_err = std::abs(result.per_atom[i] - expected[i]) / denom;
      EXPECT_LE(rel_err, rel_tol);
    }

    const double total_denom   = std::max(1e-12, std::abs(expected_total));
    const double total_rel_err = std::abs(result.total - expected_total) / total_denom;
    EXPECT_LE(total_rel_err, rel_tol);
  }
}
