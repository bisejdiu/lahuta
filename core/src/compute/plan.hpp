/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: std::invoke([] {
 *   return std::string{"besian"} + "sejdiu" + "@gmail.com";
 * });
 *
 */

#ifndef LAHUTA_COMPUTE_PLAN_HPP
#define LAHUTA_COMPUTE_PLAN_HPP

#include "node.hpp"
#include "registry.hpp"

namespace lahuta::compute {

/// Represents a validated execution plan: a sequence of computations that
/// respects dependency ordering. Validates the execution order at construction
/// time to ensure all dependencies are satisfied before their dependents.
template <typename D, Mut M = Mut::ReadWrite>
class ExecutionPlan {
  ExecOrder plan_; // immutable
public:
  explicit ExecutionPlan(const Registry<D, M> &registry, std::initializer_list<ComputationLabel> ordered) {
    // validate dependencies in O(N) time
    Mask satisfied_deps = 0;
    for (auto label : ordered) {
      int node_idx = registry.find(label);
      if (node_idx < 0) throw std::runtime_error("unknown label in execution plan");

      // every dependency of this node must already be in 'satisfied_deps'
      Mask node_deps = registry[node_idx].deps;
      for (int dep_idx = 0; dep_idx < registry.size(); ++dep_idx) {
        if (node_deps & (Mask{1} << dep_idx) && !(satisfied_deps & (Mask{1} << dep_idx))) {
          throw std::runtime_error("execution plan order violates dependency for " +
                                   std::string(registry[node_idx].tag.to_string_view()));
        }
      }

      plan_.node_indices[plan_.size++]  = static_cast<u8>(node_idx);
      satisfied_deps                   |= (Mask{1} << node_idx);
    }
  }
  const ExecOrder &plan() const noexcept { return plan_; }
};

} // namespace lahuta::compute

#endif // LAHUTA_COMPUTE_PLAN_HPP
