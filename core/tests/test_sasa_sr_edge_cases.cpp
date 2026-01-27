#include <vector>

#include <Geometry/point.h>
#include <gtest/gtest.h>

#include "analysis/sasa/sasa.hpp"
#include "utils/math_constants.hpp"
#include "utils/span.hpp"

namespace {

using lahuta::Pi;

double sphere_area(double r) { return 4.0 * Pi * r * r; }

lahuta::analysis::AtomView make_atoms(const std::vector<RDGeom::Point3D> &coords,
                                      const std::vector<double> &radii) {
  return {lahuta::span<const RDGeom::Point3D>(coords.data(), coords.size()),
          lahuta::span<const double>(radii.data(), radii.size())};
}

lahuta::analysis::AtomView make_atoms(const std::vector<RDGeom::Point3D> &coords,
                                      const std::vector<double> &radii,
                                      const std::vector<lahuta::analysis::AtomId> &ids) {
  auto view     = make_atoms(coords, radii);
  view.atom_ids = lahuta::span<const lahuta::analysis::AtomId>(ids.data(), ids.size());
  return view;
}

} // namespace

TEST(SasaSrProbeEdges, ProbeZeroMatchesArea) {
  std::vector<RDGeom::Point3D> coords = {RDGeom::Point3D(0.0, 0.0, 0.0)};
  std::vector<double> radii           = {1.5};

  lahuta::analysis::SasaParams params;
  params.probe_radius = 0.0;
  params.n_points     = 64;
  params.use_bitmask  = true;

  const auto atoms  = make_atoms(coords, radii);
  const auto result = lahuta::analysis::compute_sasa(atoms, params);

  const double expected = sphere_area(1.5);
  ASSERT_EQ(result.per_atom.size(), 1u);
  EXPECT_NEAR(result.per_atom[0], expected, 1e-12);
  EXPECT_NEAR(result.total, expected, 1e-12);
}

TEST(SasaSrProbeEdges, ProbeVeryLargeMatchesArea) {
  std::vector<RDGeom::Point3D> coords = {RDGeom::Point3D(1.0, -2.0, 3.0)};
  std::vector<double> radii           = {1.0};

  lahuta::analysis::SasaParams params;
  params.probe_radius = 100.0;
  params.n_points     = 64;
  params.use_bitmask  = true;

  const auto atoms  = make_atoms(coords, radii);
  const auto result = lahuta::analysis::compute_sasa(atoms, params);

  const double expected = sphere_area(101.0);
  ASSERT_EQ(result.per_atom.size(), 1u);
  EXPECT_NEAR(result.per_atom[0], expected, 1e-12);
  EXPECT_NEAR(result.total, expected, 1e-12);
}

TEST(SasaSrProbeEdges, TinyProbeMatchesArea) {
  std::vector<RDGeom::Point3D> coords = {RDGeom::Point3D(0.0, 0.0, 0.0)};
  std::vector<double> radii           = {1.0};

  lahuta::analysis::SasaParams params;
  params.probe_radius = 1e-6;
  params.n_points     = 64;
  params.use_bitmask  = true;

  const auto atoms  = make_atoms(coords, radii);
  const auto result = lahuta::analysis::compute_sasa(atoms, params);

  const double expected = sphere_area(1.000001);
  ASSERT_EQ(result.per_atom.size(), 1u);
  EXPECT_NEAR(result.per_atom[0], expected, 1e-12);
  EXPECT_NEAR(result.total, expected, 1e-12);
}

TEST(SasaSrProbeEdges, NegativeProbeThrows) {
  std::vector<RDGeom::Point3D> coords = {RDGeom::Point3D(0.0, 0.0, 0.0)};
  std::vector<double> radii           = {1.0};

  lahuta::analysis::SasaParams params;
  params.probe_radius = -0.1;
  params.n_points     = 64;
  params.use_bitmask  = true;

  const auto atoms = make_atoms(coords, radii);
  EXPECT_THROW({ (void)lahuta::analysis::compute_sasa(atoms, params); }, std::invalid_argument);
}

TEST(SasaSrGeometryEdges, NegativeRadiusSkipped) {
  std::vector<RDGeom::Point3D> coords = {RDGeom::Point3D(0.0, 0.0, 0.0), RDGeom::Point3D(10.0, 0.0, 0.0)};
  std::vector<double> radii           = {-1.0, 1.0};

  lahuta::analysis::SasaParams params;
  params.probe_radius = 0.0;
  params.n_points     = 64;
  params.use_bitmask  = true;

  const auto atoms          = make_atoms(coords, radii);
  const auto result_bitmask = lahuta::analysis::compute_sasa(atoms, params);

  const double expected = sphere_area(1.0);
  ASSERT_EQ(result_bitmask.per_atom.size(), 2u);
  EXPECT_NEAR(result_bitmask.per_atom[0], 0.0, 1e-12);
  EXPECT_NEAR(result_bitmask.per_atom[1], expected, 1e-12);
  EXPECT_NEAR(result_bitmask.total, expected, 1e-12);

  params.use_bitmask         = false;
  const auto result_standard = lahuta::analysis::compute_sasa(atoms, params);
  ASSERT_EQ(result_standard.per_atom.size(), 2u);
  EXPECT_NEAR(result_standard.per_atom[0], 0.0, 1e-12);
  EXPECT_NEAR(result_standard.per_atom[1], expected, 1e-12);
  EXPECT_NEAR(result_standard.total, expected, 1e-12);
}

TEST(SasaSrGeometryEdges, DuplicateAtomsSameIdNoOcclusion) {
  std::vector<RDGeom::Point3D> coords = {RDGeom::Point3D(0.0, 0.0, 0.0), RDGeom::Point3D(0.0, 0.0, 0.0)};
  std::vector<double> radii           = {1.0, 1.0};
  std::vector<lahuta::analysis::AtomId> ids = {7, 7};

  lahuta::analysis::SasaParams params;
  params.probe_radius = 0.0;
  params.n_points     = 128;
  params.use_bitmask  = true;

  const auto atoms  = make_atoms(coords, radii, ids);
  const auto result = lahuta::analysis::compute_sasa(atoms, params);

  const double expected = sphere_area(1.0);
  ASSERT_EQ(result.per_atom.size(), 2u);
  EXPECT_NEAR(result.per_atom[0], expected, 1e-12);
  EXPECT_NEAR(result.per_atom[1], expected, 1e-12);
  EXPECT_NEAR(result.total, 2.0 * expected, 1e-12);
}

TEST(SasaSrGeometryEdges, IdenticalPositionsDifferentIdsThrows) {
  std::vector<RDGeom::Point3D> coords = {RDGeom::Point3D(0.0, 0.0, 0.0), RDGeom::Point3D(0.0, 0.0, 0.0)};
  std::vector<double> radii           = {1.0, 1.0};
  std::vector<lahuta::analysis::AtomId> ids = {1, 2};

  lahuta::analysis::SasaParams params;
  params.probe_radius = 0.0;
  params.n_points     = 256;
  params.use_bitmask  = true;

  const auto atoms = make_atoms(coords, radii, ids);
  EXPECT_THROW({ (void)lahuta::analysis::compute_sasa(atoms, params); }, std::invalid_argument);

  params.use_bitmask = false;
  EXPECT_THROW({ (void)lahuta::analysis::compute_sasa(atoms, params); }, std::invalid_argument);
}

TEST(SasaSrGeometryEdges, ExtremeCoordinatesSingleSphereFinite) {
  std::vector<RDGeom::Point3D> coords = {RDGeom::Point3D(1e9, -1e9, 1e9)};
  std::vector<double> radii           = {2.5};

  lahuta::analysis::SasaParams params;
  params.probe_radius = 0.0;
  params.n_points     = 64;
  params.use_bitmask  = true;

  const auto atoms  = make_atoms(coords, radii);
  const auto result = lahuta::analysis::compute_sasa(atoms, params);

  const double expected = sphere_area(2.5);
  ASSERT_EQ(result.per_atom.size(), 1u);
  EXPECT_NEAR(result.per_atom[0], expected, 1e-12);
  EXPECT_NEAR(result.total, expected, 1e-12);
}

TEST(SasaSrGeometryEdges, MixedScaleCoordinatesFarApartFullArea) {
  std::vector<RDGeom::Point3D> coords = {RDGeom::Point3D(0.0, 0.0, 0.0), RDGeom::Point3D(1e12, 0.0, 0.0)};
  std::vector<double> radii           = {1.0, 1.0};

  lahuta::analysis::SasaParams params;
  params.probe_radius = 0.0;
  params.n_points     = 64;
  params.use_bitmask  = true;

  const auto atoms  = make_atoms(coords, radii);
  const auto result = lahuta::analysis::compute_sasa(atoms, params);

  const double expected = sphere_area(1.0);
  ASSERT_EQ(result.per_atom.size(), 2u);
  EXPECT_NEAR(result.per_atom[0], expected, 1e-12);
  EXPECT_NEAR(result.per_atom[1], expected, 1e-12);
  EXPECT_NEAR(result.total, 2.0 * expected, 1e-12);
}
