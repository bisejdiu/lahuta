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

} // namespace

TEST(SasaSrBitmask, MatchesPointsWithinTolerance) {
  std::vector<RDGeom::Point3D> coords = {
      RDGeom::Point3D(0.0, 0.0, 0.0),
      RDGeom::Point3D(3.0, 0.0, 0.0),
      RDGeom::Point3D(0.0, 3.5, 0.0),
      RDGeom::Point3D(6.0, 0.0, 0.0),
  };
  std::vector<double> radii = {2.0, 2.0, 1.7, 1.0};

  const auto atoms = make_atoms(coords, radii);

  lahuta::analysis::SasaParams params;
  params.probe_radius = 0.5;
  params.n_points     = 256;
  params.use_bitmask  = false;
  const auto baseline = lahuta::analysis::compute_sasa(atoms, params);

  params.use_bitmask = true;
  const auto bitmask = lahuta::analysis::compute_sasa(atoms, params);

  ASSERT_EQ(baseline.per_atom.size(), bitmask.per_atom.size());

  for (std::size_t i = 0; i < baseline.per_atom.size(); ++i) {
    const double full_area = sphere_area(radii[i] + params.probe_radius);
    const double tol       = full_area * 0.1;
    EXPECT_NEAR(bitmask.per_atom[i], baseline.per_atom[i], tol);
  }
  const double total_full = sphere_area(radii[0] + params.probe_radius) +
                            sphere_area(radii[1] + params.probe_radius) +
                            sphere_area(radii[2] + params.probe_radius) +
                            sphere_area(radii[3] + params.probe_radius);
  EXPECT_NEAR(bitmask.total, baseline.total, total_full * 0.1);
}
