#ifndef LAHUTA_ANALYSIS_SHAPE_METRICS_UTILS_HPP
#define LAHUTA_ANALYSIS_SHAPE_METRICS_UTILS_HPP

#include <Eigen/Dense>
#include <Geometry/point.h>

#include "utils/span.hpp"

namespace lahuta::analysis {

[[nodiscard]] Eigen::Vector3d center_of_mass(span<const RDGeom::Point3D> coords, span<const double> masses);

[[nodiscard]] Eigen::Matrix3d gyration_tensor(span<const RDGeom::Point3D> coords, span<const double> masses);

[[nodiscard]] Eigen::Vector3d principal_moments(span<const RDGeom::Point3D> coords,
                                                span<const double> masses);

[[nodiscard]] double asphericity(span<const RDGeom::Point3D> coords, span<const double> masses);

[[nodiscard]] double acylindricity(span<const RDGeom::Point3D> coords, span<const double> masses);

[[nodiscard]] double relative_shape_anisotropy(span<const RDGeom::Point3D> coords, span<const double> masses);

[[nodiscard]] double radius_of_gyration_squared(span<const RDGeom::Point3D> coords,
                                                span<const double> masses);

[[nodiscard]] double radius_of_gyration(span<const RDGeom::Point3D> coords, span<const double> masses);

} // namespace lahuta::analysis

#endif // LAHUTA_ANALYSIS_SHAPE_METRICS_UTILS_HPP
