#ifndef LAHUTA_ANALYSIS_TOPOLOGY_COMPUTATION_HPP
#define LAHUTA_ANALYSIS_TOPOLOGY_COMPUTATION_HPP

#include "analysis/system/computation.hpp"
#include "analysis/topology/kernel.hpp"
#include "compute/compute_impl.hpp"
#include "compute/dependency.hpp"
#include "pipeline/compute/context.hpp"
#include "pipeline/compute/parameters.hpp"

// clang-format off
namespace lahuta::analysis::topology {
using namespace lahuta::pipeline::compute;

class BuildTopologyComputation : public ReadWriteComputation<PipelineContext, BuildTopologyParams, BuildTopologyComputation> {
public:
  using Base = ReadWriteComputation<PipelineContext, BuildTopologyParams, BuildTopologyComputation>;
  using Base::Base;

  constexpr static const ComputationLabel label{"topology"};
  using dependencies = Dependencies<Dependency<lahuta::analysis::system::SystemReadComputation, void>>;

  ComputationResult execute_typed(DataContext<PipelineContext, Mut::ReadWrite>& ctx, const BuildTopologyParams& p) {
    return BuildTopologyKernel::execute(ctx, p);
  }
};

} // namespace lahuta::analysis::topology

#endif // LAHUTA_ANALYSIS_TOPOLOGY_COMPUTATION_HPP
