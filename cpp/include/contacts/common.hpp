#ifndef LAHUTA_CONTACTS_COMMON_HPP
#define LAHUTA_CONTACTS_COMMON_HPP

#include "Geometry/point.h"

namespace lahuta {

typedef RDGeom::Point3D Point3D;

/// Project a vector onto a plane defined by a normal vector
inline Point3D project_on_plane(const Point3D &vector, const Point3D &plane_normal) {
  // subtract component along the normal
  double scalar_projection = vector.dotProduct(plane_normal);
  RDGeom::Point3D projected_vec = vector - (plane_normal * scalar_projection);
  return projected_vec;
}

/// Compute the in-plane offset between two points
inline double compute_in_plane_offset(const Point3D &pos_a, const Point3D &pos_b, const Point3D &normal) {
  RDGeom::Point3D vec_ab = pos_a - pos_b;
  RDGeom::Point3D projected_vec = project_on_plane(vec_ab, normal);
  double in_plane_offset = projected_vec.length();
  return in_plane_offset;
}

} // namespace lahuta

#endif // LAHUTA_COMMON_HPP
