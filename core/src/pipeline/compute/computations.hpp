#ifndef LAHUTA_PIPELINE_COMPUTE_COMPUTATIONS_HPP
#define LAHUTA_PIPELINE_COMPUTE_COMPUTATIONS_HPP

#include <memory>
#include <string>
#include <vector>

#include "compute/compute_base.hpp"
#include "compute/result.hpp"
#include "pipeline/compute/context.hpp"
#include "pipeline/compute/parameters.hpp"
#include "pipeline/dynamic/types.hpp"

// clang-format off
namespace lahuta::pipeline::compute {
using namespace lahuta::topology::compute;

// Dynamic RW base used for dynamic tasks: runtime label and runtime deps
template <typename P>
class DynamicRWComputation : public Computation<PipelineContext, Mut::ReadWrite> {
public:
  DynamicRWComputation(std::string label, std::vector<std::string> dep_names, P params)
      : label_storage_(std::move(label)), label_(label_storage_), dep_names_(std::move(dep_names)), params_(std::move(params)) {}

  ComputationResult execute(DataContext<PipelineContext, Mut::ReadWrite>& context, const ParameterInterface& raw) override {
    if (raw.type_id() != P::TYPE_ID) {
      return ComputationResult(ComputationError("Invalid parameter type for runtime computation: " + std::string(label_.to_string_view())));
    }
    const auto& typed = static_cast<const P&>(raw);
    return execute_typed(context, typed);
  }

  std::unique_ptr<ParameterInterface> get_parameters() const override { return std::make_unique<P>(params_); }
  const ComputationLabel& get_label() const override { return label_; }
  std::vector<ComputationLabel> get_dependencies() const override {
    std::vector<ComputationLabel> deps;
    deps.reserve(dep_names_.size());
    for (const auto& s : dep_names_) deps.emplace_back(s.c_str());
    return deps;
  }

protected:
  virtual ComputationResult execute_typed(DataContext<PipelineContext, Mut::ReadWrite>&, const P&) = 0;

private:
  std::string label_storage_;
  ComputationLabel label_;
  std::vector<std::string> dep_names_;
  P params_;
};

//
// Bridge for dynamic tasks
//
class DynamicTaskComputation final : public DynamicRWComputation<DynamicTaskParams> {
public:
  DynamicTaskComputation(std::string label, std::vector<std::string> dep_names, std::shared_ptr<dynamic::ITask> task)
    : DynamicRWComputation<DynamicTaskParams>(std::move(label), std::move(dep_names), DynamicTaskParams{}), task_(std::move(task)) {}

  pipeline::DataFieldSet data_requirements() const override {
    return task_ ? task_->data_requirements() : pipeline::DataFieldSet::none();
  }

private:
  ComputationResult execute_typed(DataContext<PipelineContext, Mut::ReadWrite>& context, const DynamicTaskParams&) override {
    auto& data = context.data();
    if (!task_) return ComputationResult(ComputationError("DynamicTaskComputation: null task"));

    auto r = task_->run(data.item_path, *data.ctx);
    if (!r.ok) return ComputationResult(ComputationError("Dynamic task returned ok=false"));

    dynamic::EmissionList out;
    out.reserve(r.emits.size());
    for (auto &e : r.emits) out.push_back(std::move(e));
    return ComputationResult(std::move(out));
  }

  std::shared_ptr<dynamic::ITask> task_;
};

} // namespace lahuta::pipeline::compute

#endif // LAHUTA_PIPELINE_COMPUTE_COMPUTATIONS_HPP
