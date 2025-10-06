#ifndef LAHUTA_DISTANCES_NEIGHBORS_HPP
#define LAHUTA_DISTANCES_NEIGHBORS_HPP

#include <vector>

#include "nsgrid.hpp"
#include <rdkit/Geometry/point.h>

// clang-format off
namespace lahuta::dist {

struct NeighborSearchOptions {
  double cutoff{0.0};
  bool brute_force_fallback{true};
  bool sort_output{false};
};

NSResults neighbors_within_radius_self(const RDGeom::POINT3D_VECT &coords, const NeighborSearchOptions &options);
NSResults neighbors_within_radius_cross(const RDGeom::POINT3D_VECT &queries, const RDGeom::POINT3D_VECT &targets, const NeighborSearchOptions &options);
NSResults neighbors_within_radius_cross_fastns(const RDGeom::POINT3D_VECT &queries, const RDGeom::POINT3D_VECT &targets, const NeighborSearchOptions &options);

NSResults brute_force_radius_self_streamed(const RDGeom::POINT3D_VECT &coords, double cutoff);

NSResults brute_force_radius_cross_streamed(const RDGeom::POINT3D_VECT &queries, const RDGeom::POINT3D_VECT &targets, double cutoff);

template <typename Real>
std::vector<Real> to_euclidean_distances(const NSResults &results) {
  std::vector<Real> out;
  out.reserve(results.get_distances().size());
  for (float d2 : results.get_distances()) {
    out.push_back(static_cast<Real>(std::sqrt(static_cast<Real>(d2))));
  }
  return out;
}

} // namespace lahuta::dist

#endif // LAHUTA_DISTANCES_NEIGHBORS_HPP
