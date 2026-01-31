#ifndef LAHUTA_ANALYSIS_SASA_METHODS_BITMASK_HPP
#define LAHUTA_ANALYSIS_SASA_METHODS_BITMASK_HPP

#include "analysis/sasa/bitmask_lut.hpp"
#include "analysis/sasa/method.hpp"

namespace lahuta::analysis {

// Bitmask-based SASA method using precomputed lookup tables.
// Uses bitwise operations for occlusion testing, much faster than point-by-point.
class BitmaskMethod final : public SasaMethod {
public:
  explicit BitmaskMethod(const BitmaskLut &lut) : lut_(lut) {}

  [[nodiscard]] double compute(const ComputeContext<CoordAccessor> &ctx) const override;

private:
  template <typename CoordAccessor>
  [[nodiscard]] double compute_impl(const ComputeContext<CoordAccessor> &ctx) const;

  const BitmaskLut &lut_;
};

} // namespace lahuta::analysis

#endif // LAHUTA_ANALYSIS_SASA_METHODS_BITMASK_HPP
