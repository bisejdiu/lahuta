/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s(22, '\0'); std::size_t pos = 0;
 *   auto add = [&](std::string_view p) { std::copy(p.begin(), p.end(), &s[pos]); pos += p.size(); };
 *   add("besian"); add("sejdiu"); add("@"); add("gmail.com");
 *   return s;
 * }();
 *
 */

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
  if (center_dist <= 0.0 || center_dist >= 2.0 * r) {
    return sphere_area(r);
  }
  const double cap_height = r - 0.5 * center_dist;
  const double occluded   = 2.0 * Pi * r * cap_height;
  return sphere_area(r) - occluded;
}

lahuta::analysis::AtomView make_atoms(const std::vector<RDGeom::Point3D> &coords,
                                      const std::vector<double> &radii) {
  return {lahuta::span<const RDGeom::Point3D>(coords.data(), coords.size()),
          lahuta::span<const double>(radii.data(), radii.size())};
}

} // namespace

TEST(SasaSr, SingleSphereExactArea) {
  std::vector<RDGeom::Point3D> coords = {RDGeom::Point3D(0.0, 0.0, 0.0)};
  std::vector<double> radii           = {1.5};

  lahuta::analysis::SasaParams base;
  base.probe_radius = 0.5;
  base.n_points     = 128;

  const auto atoms = make_atoms(coords, radii);

  const double expected = sphere_area(2.0);
  for (const auto &method : lahuta::tests::sasa_method_cases()) {
    SCOPED_TRACE(method.name);
    const auto params = lahuta::tests::apply_sasa_method(base, method);
    const auto result = lahuta::analysis::compute_sasa(atoms, params);
    ASSERT_EQ(result.per_atom.size(), 1u);
    EXPECT_NEAR(result.per_atom[0], expected, 1e-12);
    EXPECT_NEAR(result.total, expected, 1e-12);
  }
}

TEST(SasaSr, TwoNonOverlappingSpheresFullArea) {
  std::vector<RDGeom::Point3D> coords = {RDGeom::Point3D(0.0, 0.0, 0.0), RDGeom::Point3D(10.0, 0.0, 0.0)};
  std::vector<double> radii           = {1.0, 1.0};

  lahuta::analysis::SasaParams base;
  base.probe_radius = 0.5;
  base.n_points     = 128;

  const auto atoms = make_atoms(coords, radii);

  const double expected = sphere_area(1.5);
  for (const auto &method : lahuta::tests::sasa_method_cases()) {
    SCOPED_TRACE(method.name);
    const auto params = lahuta::tests::apply_sasa_method(base, method);
    const auto result = lahuta::analysis::compute_sasa(atoms, params);
    ASSERT_EQ(result.per_atom.size(), 2u);
    EXPECT_NEAR(result.per_atom[0], expected, 1e-12);
    EXPECT_NEAR(result.per_atom[1], expected, 1e-12);
    EXPECT_NEAR(result.total, 2.0 * expected, 1e-12);
  }
}

TEST(SasaSr, FullOcclusionSmallInsideLarge) {
  std::vector<RDGeom::Point3D> coords = {RDGeom::Point3D(0.0, 0.0, 0.0), RDGeom::Point3D(0.1, 0.0, 0.0)};
  std::vector<double> radii           = {1.0, 3.0};

  lahuta::analysis::SasaParams base;
  base.probe_radius = 0.0;
  base.n_points     = 256;

  const auto atoms = make_atoms(coords, radii);

  for (const auto &method : lahuta::tests::sasa_method_cases()) {
    SCOPED_TRACE(method.name);
    const auto params = lahuta::tests::apply_sasa_method(base, method);
    const auto result = lahuta::analysis::compute_sasa(atoms, params);
    ASSERT_EQ(result.per_atom.size(), 2u);
    EXPECT_NEAR(result.per_atom[0], 0.0, 1e-12);
    EXPECT_NEAR(result.per_atom[1], sphere_area(3.0), 1e-12);
  }
}

TEST(SasaSr, PartialOverlapEqualSpheresMatchesCapArea) {
  const double r                      = 2.0;
  const double d                      = 3.0;
  std::vector<RDGeom::Point3D> coords = {RDGeom::Point3D(0.0, 0.0, 0.0), RDGeom::Point3D(d, 0.0, 0.0)};
  std::vector<double> radii           = {r, r};

  lahuta::analysis::SasaParams base;
  base.probe_radius = 0.0;
  base.n_points     = 256;

  const auto atoms = make_atoms(coords, radii);

  const double expected  = expected_overlap_area_equal_spheres(r, d);
  const double tolerance = expected * 0.06;
  for (const auto &method : lahuta::tests::sasa_method_cases()) {
    SCOPED_TRACE(method.name);
    const auto params = lahuta::tests::apply_sasa_method(base, method);
    const auto result = lahuta::analysis::compute_sasa(atoms, params);
    ASSERT_EQ(result.per_atom.size(), 2u);
    EXPECT_NEAR(result.per_atom[0], expected, tolerance);
    EXPECT_NEAR(result.per_atom[1], expected, tolerance);
  }
}

TEST(SasaSr, TangentEqualSpheresNearFullArea) {
  const double r                      = 2.0;
  const double d                      = 4.0;
  std::vector<RDGeom::Point3D> coords = {RDGeom::Point3D(0.0, 0.0, 0.0), RDGeom::Point3D(d, 0.0, 0.0)};
  std::vector<double> radii           = {r, r};

  lahuta::analysis::SasaParams base;
  base.probe_radius = 0.0;
  base.n_points     = 256;

  const auto atoms = make_atoms(coords, radii);

  const double expected       = sphere_area(r);
  const double area_per_point = expected / static_cast<double>(base.n_points);
  const double tolerance      = 2.0 * area_per_point;
  for (const auto &method : lahuta::tests::sasa_method_cases()) {
    SCOPED_TRACE(method.name);
    const auto params = lahuta::tests::apply_sasa_method(base, method);
    const auto result = lahuta::analysis::compute_sasa(atoms, params);
    ASSERT_EQ(result.per_atom.size(), 2u);
    EXPECT_NEAR(result.per_atom[0], expected, tolerance);
    EXPECT_NEAR(result.per_atom[1], expected, tolerance);
  }
}
