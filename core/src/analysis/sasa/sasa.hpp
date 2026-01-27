#ifndef LAHUTA_ANALYSIS_SASA_SASA_HPP
#define LAHUTA_ANALYSIS_SASA_SASA_HPP

#include <cstddef>
#include <vector>

#include <Geometry/point.h>

#include "analysis/sasa/types.hpp"
#include "utils/span.hpp"

namespace lahuta::analysis {

struct SasaParams {
  double probe_radius       = 1.4;
  SpherePointCount n_points = 100;
  bool skip_same_id         = true;
  bool use_bitmask          = true;
};

struct AtomView {
  span<const RDGeom::Point3D> coords;
  span<const double> radii;
  span<const AtomId> atom_ids{};
};

struct SasaResult {
  std::vector<double> per_atom;    // SASA per atom in A^2
  double total              = 0.0; // Total SASA in A^2
  double probe_radius       = 0.0;
  SpherePointCount n_points = 0;
};

// Compute SASA using Shrake-Rupley
[[nodiscard]] SasaResult compute_sasa(const AtomView &atoms, const SasaParams &params);

} // namespace lahuta::analysis

#endif // LAHUTA_ANALYSIS_SASA_SASA_HPP
