#ifndef LAHUTA_PIPELINE_COMPUTE_DYNAMIC_COMPUTATION_HPP
#define LAHUTA_PIPELINE_COMPUTE_DYNAMIC_COMPUTATION_HPP

#include <string>
#include <vector>

#include "compute/compute_base.hpp"
#include "pipeline/compute/context.hpp"

// clang-format off
namespace lahuta::pipeline::compute {
using namespace lahuta::topology::compute;

// Computations that need runtime labels but have compile-time dependencies.
template <typename P, typename KernelT, typename Derived>
class DynamicLabelComputation : public Computation<PipelineContext, Mut::ReadWrite> {
public:
  explicit DynamicLabelComputation(std::string label, P params)
      : label_storage_(std::move(label)), label_(label_storage_), params_(std::move(params)) {}

  ComputationResult execute(DataContext<PipelineContext, Mut::ReadWrite>& ctx, const ParameterInterface& raw) override {
    if (raw.type_id() != P::TYPE_ID) {
      return ComputationResult(ComputationError("Invalid parameter type for " + std::string(label_.to_string_view())));
    }
    const auto& typed = static_cast<const P&>(raw);
    return KernelT::execute(ctx, typed);
  }

  std::unique_ptr<ParameterInterface> get_parameters() const override { return std::make_unique<P>(params_); }
  const ComputationLabel& get_label() const override { return label_; }
  std::vector<ComputationLabel> get_dependencies() const override {
    return compute::detail::dependencies_of<Derived>::labels();
  }

private:
  std::string label_storage_;
  ComputationLabel label_;
  P params_;
};

} // namespace lahuta::pipeline::compute

#endif // LAHUTA_PIPELINE_COMPUTE_DYNAMIC_COMPUTATION_HPP
