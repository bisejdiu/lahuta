#pragma once

#include "compute/_defs.hpp"
#include "node.hpp"
#include "registry.hpp"

// clang-format off
namespace lahuta::topology::compute {

/// Schedules and executes computations with dependency tracking. Returns the execution order of computations.
// - determines execution order
// - computes, respects enabled/disabled identifiers, memoizes results
// - returns the execution order, throws an error if a cycle is detected
template<typename D, Mut M>
ExecOrder schedule_and_run(Registry<D, M>& registry, DataContext<D, M>& ctx) {
  ExecOrder exec_order{};
  Mask satisfied_deps = 0, in_progress = 0;

  auto dfs = [&](auto& self, int node_idx) -> void {
    if(!registry[node_idx].enabled) return;
    if(satisfied_deps    & (1<<node_idx)) return;
    if(in_progress       & (1<<node_idx)) throw std::runtime_error("cycle @ " + std::string(registry[node_idx].tag.to_string_view()));

    in_progress |= (1<<node_idx);
    for(int dep_idx=0; dep_idx<registry.size(); ++dep_idx) {
      if(registry[node_idx].deps & (1<<dep_idx)) {
        self(self, dep_idx);
      }
    }
    in_progress &= ~(1<<node_idx);

    // all deps executed now
    auto& node = registry[node_idx];
    if(!node.done) {
      auto params = node.proto->clone();
      node.res    = node.impl->execute(ctx, *params);
      node.done   = true;

    }
    exec_order.node_indices[exec_order.size++] = node_idx;
    satisfied_deps |= (1<<node_idx);
  };

  for(int node_idx=0; node_idx<registry.size(); ++node_idx) {
    if(registry[node_idx].enabled) dfs(dfs, node_idx);
  }
  return exec_order;
}

/// Executes computations given by the pipeline. It is the caller's responsibility to ensure that the pipeline is valid.
// - assumes the pipeline is already validated and sealed
// - does not check for cycles or dependencies
// - does not check for enabled/disabled identifiers.
// - does memoize results
template<typename D, Mut M>
void execute_pipeline(Registry<D, M>& registry, DataContext<D, M>& ctx, const ExecOrder& exec_plan) {
  for (u8 plan_idx = 0; plan_idx < exec_plan.size; ++plan_idx) {
    ComputeNode<D, M>& node = registry[exec_plan.node_indices[plan_idx]];
    if (node.done) continue; // memoized

    auto params = node.proto->clone();
    node.res    = node.impl->execute(ctx, *params);
    node.done   = true;
  }
}

} // namespace lahuta::topology::compute
