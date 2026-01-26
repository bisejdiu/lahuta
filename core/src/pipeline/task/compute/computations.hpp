#ifndef LAHUTA_PIPELINE_TASK_COMPUTE_COMPUTATIONS_HPP
#define LAHUTA_PIPELINE_TASK_COMPUTE_COMPUTATIONS_HPP

#include <memory>
#include <string>
#include <vector>

#include "compute/compute_base.hpp"
#include "compute/result.hpp"
#include "pipeline/task/compute/context.hpp"
#include "pipeline/task/compute/parameters.hpp"
#include "pipeline/task/task.hpp"

namespace lahuta::pipeline {
namespace C = lahuta::compute;
namespace P = lahuta::pipeline;

// Dynamic RW base used for dynamic tasks: runtime label and runtime deps
template <typename P>
class DynamicRWComputation : public C::Computation<PipelineContext> {
public:
  DynamicRWComputation(std::string label, std::vector<std::string> dep_names, P params)
      : label_storage_(std::move(label)), label_(label_storage_), dep_names_(std::move(dep_names)),
        params_(std::move(params)) {}

  C::ComputationResult execute(C::DataContext<PipelineContext> &context,
                               const C::ParameterInterface &raw) override {
    if (raw.type_id() != P::TYPE_ID) {
      return C::ComputationResult(C::ComputationError("Invalid parameter type for runtime computation: " +
                                                      std::string(label_.to_string_view())));
    }
    const auto &typed = static_cast<const P &>(raw);
    return execute_typed(context, typed);
  }

  std::unique_ptr<C::ParameterInterface> get_parameters() const override {
    return std::make_unique<P>(params_);
  }
  const C::ComputationLabel &get_label() const override { return label_; }
  std::vector<C::ComputationLabel> get_dependencies() const override {
    std::vector<C::ComputationLabel> deps;
    deps.reserve(dep_names_.size());
    for (const auto &s : dep_names_)
      deps.emplace_back(s.c_str());
    return deps;
  }

protected:
  virtual C::ComputationResult execute_typed(C::DataContext<PipelineContext> &, const P &) = 0;

private:
  std::string label_storage_;
  C::ComputationLabel label_;
  std::vector<std::string> dep_names_;
  P params_;
};

//
// Bridge for dynamic tasks
//
class DynamicTaskComputation final : public DynamicRWComputation<DynamicTaskParams> {
public:
  DynamicTaskComputation(std::string label, std::vector<std::string> dep_names,
                         std::shared_ptr<P::ITask> task)
      : DynamicRWComputation<DynamicTaskParams>(std::move(label), std::move(dep_names), DynamicTaskParams{}),
        task_(std::move(task)) {}

  P::DataFieldSet data_requirements() const override {
    return task_ ? task_->data_requirements() : P::DataFieldSet::none();
  }

private:
  C::ComputationResult execute_typed(C::DataContext<PipelineContext> &context,
                                     const DynamicTaskParams &) override {
    auto &data = context.data();
    if (!task_) return C::ComputationResult(C::ComputationError("DynamicTaskComputation: null task"));

    auto r = task_->run(data.item_path, *data.ctx);
    if (!r.ok) return C::ComputationResult(C::ComputationError("Dynamic task returned ok=false"));

    P::EmissionList out;
    out.reserve(r.emits.size());
    for (auto &e : r.emits)
      out.push_back(std::move(e));
    return C::ComputationResult(std::move(out));
  }

  std::shared_ptr<P::ITask> task_;
};

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_TASK_COMPUTE_COMPUTATIONS_HPP
