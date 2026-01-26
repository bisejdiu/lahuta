#ifndef LAHUTA_ANALYSIS_TOPOLOGY_COMPUTATION_HPP
#define LAHUTA_ANALYSIS_TOPOLOGY_COMPUTATION_HPP

#include "analysis/system/computation.hpp"
#include "analysis/topology/kernel.hpp"
#include "compute/compute_impl.hpp"
#include "compute/dependency.hpp"
#include "pipeline/task/compute/context.hpp"
#include "pipeline/task/compute/parameters.hpp"

namespace lahuta::analysis {
namespace C = lahuta::compute;
namespace P = lahuta::pipeline;

class BuildTopologyComputation
    : public C::ReadWriteComputation<P::PipelineContext, P::BuildTopologyParams, BuildTopologyComputation> {
public:
  using Base = C::ReadWriteComputation<P::PipelineContext, P::BuildTopologyParams, BuildTopologyComputation>;
  using Base::Base;

  constexpr static const C::ComputationLabel label{"topology"};
  using dependencies = C::Dependencies<C::Dependency<SystemReadComputation, void>>;

  using RWContext = C::DataContext<P::PipelineContext, C::Mut::ReadWrite>;

  C::ComputationResult execute_typed(RWContext &ctx, const P::BuildTopologyParams &p) {
    return BuildTopologyKernel::execute(ctx, p);
  }
};

} // namespace lahuta::analysis

#endif // LAHUTA_ANALYSIS_TOPOLOGY_COMPUTATION_HPP
