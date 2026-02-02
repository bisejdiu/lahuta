#ifndef LAHUTA_ANALYSIS_SASA_COMPUTATION_HPP
#define LAHUTA_ANALYSIS_SASA_COMPUTATION_HPP

#include "analysis/sasa/kernel.hpp"
#include "pipeline/task/compute/dynamic_computation.hpp"

namespace lahuta::analysis {
namespace P = lahuta::pipeline;

class SasaSrComputation
    : public P::DynamicLabelComputation<P::SasaSrParams, SasaSrKernel, SasaSrComputation> {
public:
  using Base = P::DynamicLabelComputation<P::SasaSrParams, SasaSrKernel, SasaSrComputation>;
  using Base::DynamicLabelComputation;

  P::DataFieldSet data_requirements() const override {
    return P::DataFieldSet::of({P::DataField::Sequence, P::DataField::Positions});
  }
};

} // namespace lahuta::analysis

#endif // LAHUTA_ANALYSIS_SASA_COMPUTATION_HPP
