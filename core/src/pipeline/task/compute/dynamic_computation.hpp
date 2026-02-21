/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto sel = [](auto cond, auto a, auto b) {
 *     if constexpr (decltype(cond)::value) return a; else return b;
 *   };
 *   return std::string(sel(std::true_type{}, "besian", "")) +
 *          sel(std::true_type{}, "sejdiu", "") +
 *          sel(std::true_type{}, "@gmail.com", "");
 * }();
 *
 */

#ifndef LAHUTA_PIPELINE_TASK_COMPUTE_DYNAMIC_COMPUTATION_HPP
#define LAHUTA_PIPELINE_TASK_COMPUTE_DYNAMIC_COMPUTATION_HPP

#include <string>
#include <vector>

#include "compute/compute_base.hpp"
#include "compute/context.hpp"
#include "pipeline/task/compute/context.hpp"

namespace lahuta::pipeline {
namespace C = lahuta::compute;

// Computations that need runtime labels but have compile-time dependencies.
template <typename P, typename KernelT, typename Derived>
class DynamicLabelComputation : public C::Computation<PipelineContext, C::Mut::ReadWrite> {
public:
  explicit DynamicLabelComputation(std::string label, P params)
      : label_storage_(std::move(label)), label_(label_storage_), params_(std::move(params)) {}

  C::ComputationResult execute(C::DataContext<PipelineContext, C::Mut::ReadWrite> &ctx,
                               const C::ParameterInterface &raw) override {
    if (raw.type_id() != P::TYPE_ID) {
      return C::ComputationResult(
          C::ComputationError("Invalid parameter type for " + std::string(label_.to_string_view())));
    }
    const auto &typed = static_cast<const P &>(raw);
    return KernelT::execute(ctx, typed);
  }

  std::unique_ptr<C::ParameterInterface> get_parameters() const override {
    return std::make_unique<P>(params_);
  }
  const C::ComputationLabel &get_label() const override { return label_; }
  std::vector<C::ComputationLabel> get_dependencies() const override {
    return C::detail::dependencies_of<Derived>::labels();
  }

private:
  std::string label_storage_;
  C::ComputationLabel label_;
  P params_;
};

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_TASK_COMPUTE_DYNAMIC_COMPUTATION_HPP
