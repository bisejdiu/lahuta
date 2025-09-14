#ifndef LAHUTA_ARPEGGIO_GEO_HPP
#define LAHUTA_ARPEGGIO_GEO_HPP

#include <entities/records.hpp>

// clang-format off
namespace lahuta::arpeggio {

inline double compute_angle(const RingRec &rd, const RDGeom::Point3D &point) {
  auto vector_point_to_plane = point - rd.center;
  vector_point_to_plane.normalize();

  double cos_theta = vector_point_to_plane.dotProduct(rd.normal);
  cos_theta = std::max(-1.0, std::min(1.0, cos_theta));

  double theta_radians = std::acos(cos_theta); // in radians
  return theta_radians * (180.0 / M_PI);
}

inline bool passes_angle_filter(double angle, double cutoff) {
  if (angle > 90.0) angle = 180.0 - angle; // fold angles >90 back in [0,90]
  return angle <= cutoff;
}

} // namespace lahuta::arpeggio

#endif // LAHUTA_ARPEGGIO_GEO_HPP
