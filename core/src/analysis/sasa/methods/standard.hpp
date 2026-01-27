#ifndef LAHUTA_ANALYSIS_SASA_METHODS_STANDARD_HPP
#define LAHUTA_ANALYSIS_SASA_METHODS_STANDARD_HPP

#include "analysis/sasa/method.hpp"
#include "analysis/sasa/sphere.hpp"

namespace lahuta::analysis {

// Standard Shrake-Rupley SASA method using point-by-point visibility testing.
class StandardMethod final : public SasaMethod {
public:
  explicit StandardMethod(const SpherePoints &sphere) : sphere_(sphere) {}

  [[nodiscard]] double compute(const ComputeContext<CoordAccessor> &ctx) const override;

private:
  template <typename CoordAccessor>
  [[nodiscard]] double compute_impl(const ComputeContext<CoordAccessor> &ctx) const;

  const SpherePoints &sphere_;
};

} // namespace lahuta::analysis

#endif // LAHUTA_ANALYSIS_SASA_METHODS_STANDARD_HPP
