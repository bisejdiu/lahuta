#ifndef LAHUTA_ANALYSIS_SASA_METHOD_HPP
#define LAHUTA_ANALYSIS_SASA_METHOD_HPP

#include <cassert>
#include <cstddef>

#include "analysis/sasa/neighbor.hpp"
#include "utils/span.hpp"

namespace lahuta::analysis {

// Context passed to SASA methods for per-atom computation.
template <typename CoordAccessor>
struct ComputeContext {
  AtomIndex atom_index;
  const CoordAccessor &coords;
  const NeighborList &neighbors;
  span<const double> expanded_radii;
  span<const double> expanded_radii_sq;

  [[nodiscard]] std::size_t neighbor_start() const noexcept {
    assert(atom_index < neighbors.offsets.size());
    return neighbors.offsets[atom_index];
  }

  [[nodiscard]] std::size_t neighbor_end() const noexcept {
    assert(atom_index + 1 < neighbors.offsets.size());
    return neighbors.offsets[atom_index + 1];
  }

  [[nodiscard]] double radius() const noexcept { return expanded_radii[atom_index]; }
  [[nodiscard]] double radius_sq() const noexcept { return expanded_radii_sq[atom_index]; }
};

// Methods:
// - StandardMethod: Point-by-point visibility test (Shrake-Rupley)
// - BitmaskMethod: LUT-based bitmask occlusion
// - Local Overlap Polynomial SASA: attempted, but not enabled
//
// I tried implementing LOP-SASA, and while it was ~30-50% faster than the bitmask
// method, accuracy was noticeably worse. May revisit this in the future if there
// is interest, but unclear if there's any advantage over simple MLP.  - Besian, Jan 2026
class SasaMethod {
public:
  virtual ~SasaMethod() = default;

  [[nodiscard]] virtual double compute(const ComputeContext<CoordAccessor> &ctx) const = 0;
};

} // namespace lahuta::analysis

#endif // LAHUTA_ANALYSIS_SASA_METHOD_HPP
