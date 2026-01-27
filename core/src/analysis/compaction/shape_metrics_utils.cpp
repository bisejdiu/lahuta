#include <cmath>
#include <cstddef>
#include <stdexcept>

#include "analysis/compaction/shape_metrics_utils.hpp"

namespace lahuta::analysis {
namespace {

struct MassStats {
  Eigen::Vector3d center = Eigen::Vector3d::Zero();
  double total_mass      = 0.0;
};

MassStats compute_mass_stats(span<const RDGeom::Point3D> coords, span<const double> masses) {
  if (coords.size() != masses.size()) {
    throw std::runtime_error("Expected masses to match coordinate count in shape metrics");
  }
  MassStats stats{};
  if (coords.empty()) return stats;

  for (std::size_t i = 0; i < coords.size(); ++i) {
    const double mass  = masses[i];
    stats.total_mass  += mass;
    stats.center      += mass * Eigen::Vector3d{coords[i].x, coords[i].y, coords[i].z};
  }

  if (stats.total_mass <= 0.0) {
    stats.center.setZero();
    stats.total_mass = 0.0;
    return stats;
  }

  stats.center /= stats.total_mass;
  return stats;
}

} // namespace

Eigen::Vector3d center_of_mass(span<const RDGeom::Point3D> coords, span<const double> masses) {
  return compute_mass_stats(coords, masses).center;
}

Eigen::Matrix3d gyration_tensor(span<const RDGeom::Point3D> coords, span<const double> masses) {
  if (coords.empty()) return Eigen::Matrix3d::Zero();

  const auto stats = compute_mass_stats(coords, masses);
  if (stats.total_mass <= 0.0) return Eigen::Matrix3d::Zero();

  Eigen::Matrix3d tensor = Eigen::Matrix3d::Zero();
  for (std::size_t i = 0; i < coords.size(); ++i) {
    const Eigen::Vector3d r{coords[i].x - stats.center.x(),
                            coords[i].y - stats.center.y(),
                            coords[i].z - stats.center.z()};
    tensor += masses[i] * (r * r.transpose());
  }

  return tensor / stats.total_mass;
}

Eigen::Vector3d principal_moments(span<const RDGeom::Point3D> coords, span<const double> masses) {
  const Eigen::Matrix3d tensor = gyration_tensor(coords, masses);
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(tensor, Eigen::EigenvaluesOnly);
  if (solver.info() != Eigen::Success) {
    throw std::runtime_error("Failed to compute principal moments from gyration tensor");
  }
  return solver.eigenvalues();
}

double asphericity(span<const RDGeom::Point3D> coords, span<const double> masses) {
  const Eigen::Vector3d pm = principal_moments(coords, masses);
  return pm[2] - 0.5 * (pm[0] + pm[1]);
}

double acylindricity(span<const RDGeom::Point3D> coords, span<const double> masses) {
  const Eigen::Vector3d pm = principal_moments(coords, masses);
  return pm[1] - pm[0];
}

double relative_shape_anisotropy(span<const RDGeom::Point3D> coords, span<const double> masses) {
  const Eigen::Vector3d pm = principal_moments(coords, masses);
  const double sum         = pm.sum();
  if (sum <= 0.0) return 0.0;
  return 1.5 * pm.squaredNorm() / (sum * sum) - 0.5;
}

double radius_of_gyration_squared(span<const RDGeom::Point3D> coords, span<const double> masses) {
  return principal_moments(coords, masses).sum();
}

double radius_of_gyration(span<const RDGeom::Point3D> coords, span<const double> masses) {
  const double rg_sq = radius_of_gyration_squared(coords, masses);
  return rg_sq > 0.0 ? std::sqrt(rg_sq) : 0.0;
}

} // namespace lahuta::analysis
