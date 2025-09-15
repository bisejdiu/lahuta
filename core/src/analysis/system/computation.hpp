#pragma once

#include "analysis/system/kernel.hpp"
#include "compute/compute_impl.hpp"
#include "compute/dependency.hpp"
#include "pipeline/compute/context.hpp"
#include "pipeline/compute/parameters.hpp"
#include "pipeline/compute/dynamic_computation.hpp"
#include <analysis/system/model_fetch.hpp>

// clang-format off
namespace lahuta::analysis::system {
using namespace lahuta::topology::compute;

class SystemReadComputation : public ReadWriteComputation<PipelineContext, SystemReadParams, SystemReadComputation> {
public:
  using Base = ReadWriteComputation<PipelineContext, SystemReadParams, SystemReadComputation>;
  using Base::Base;

  constexpr static const ComputationLabel label{"system"};
  using dependencies = UnitComputation;

  ComputationResult execute_typed(DataContext<PipelineContext, Mut::ReadWrite>& ctx, const SystemReadParams& p) { 
    return SystemReadKernel::execute(ctx, p);
  }
};

class ModelFetchComputation : public DynamicLabelComputation<ModelFetchParams, ModelFetchKernel, ModelFetchComputation> {
public:
  using Base = DynamicLabelComputation<ModelFetchParams, ModelFetchKernel, ModelFetchComputation>;
  using Base::DynamicLabelComputation;
  using dependencies = UnitComputation;
};

} // namespace lahuta::analysis::system
