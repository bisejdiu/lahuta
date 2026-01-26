#ifndef LAHUTA_ANALYSIS_SYSTEM_COMPUTATION_HPP
#define LAHUTA_ANALYSIS_SYSTEM_COMPUTATION_HPP

#include "analysis/system/kernel.hpp"
#include "compute/compute_impl.hpp"
#include "compute/dependency.hpp"
#include "pipeline/data/data_requirements.hpp"
#include "pipeline/task/compute/context.hpp"
#include "pipeline/task/compute/parameters.hpp"

namespace lahuta::analysis {
namespace C = lahuta::compute;
namespace P = lahuta::pipeline;

class SystemReadComputation
    : public C::ReadWriteComputation<P::PipelineContext, P::SystemReadParams, SystemReadComputation> {
public:
  using Base = C::ReadWriteComputation<P::PipelineContext, P::SystemReadParams, SystemReadComputation>;
  using Base::Base;

  constexpr static const C::ComputationLabel label{"system"};
  using dependencies = C::UnitComputation;

  using RWContext = C::DataContext<P::PipelineContext, C::Mut::ReadWrite>;

  C::ComputationResult execute_typed(RWContext &ctx, const P::SystemReadParams &p) {
    return SystemReadKernel::execute(ctx, p);
  }

  P::DataFieldSet data_requirements() const override {
    return P::DataFieldSet::of({P::DataField::Sequence, P::DataField::Positions, P::DataField::Plddt});
  }
};

} // namespace lahuta::analysis

#endif // LAHUTA_ANALYSIS_SYSTEM_COMPUTATION_HPP
