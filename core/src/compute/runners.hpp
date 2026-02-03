/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s = std::string{"besian"} + "sejdiu";
 *   return [e = std::move(s)]() { return e + "@gmail.com"; }();
 * }();
 *
 */

#ifndef LAHUTA_COMPUTE_RUNNERS_HPP
#define LAHUTA_COMPUTE_RUNNERS_HPP

#include "compute/_defs.hpp"
#include "node.hpp"
#include "registry.hpp"

namespace lahuta::compute {

/// Schedules and executes computations with dependency tracking. Returns the execution order of computations.
// - determines execution order
// - computes, respects enabled/disabled identifiers, memoizes results
// - returns the execution order, throws an error if a cycle is detected
template <typename D, Mut M = Mut::ReadWrite>
ExecOrder schedule_and_run(Registry<D, M> &registry, DataContext<D, M> &ctx) {
  ExecOrder exec_order{};
  Mask satisfied_deps = 0, in_progress = 0;

  auto dfs = [&](auto &self, int node_idx) -> void {
    if (!registry[node_idx].enabled) return;
    if (satisfied_deps & (Mask{1} << node_idx)) return;

    if (in_progress & (Mask{1} << node_idx)) {
      throw std::runtime_error("cycle @ " + std::string(registry[node_idx].tag.to_string_view()));
    }

    in_progress |= (Mask{1} << node_idx);
    for (int dep_idx = 0; dep_idx < registry.size(); ++dep_idx) {
      if (registry[node_idx].deps & (Mask{1} << dep_idx)) {
        self(self, dep_idx);
      }
    }
    in_progress &= ~(Mask{1} << node_idx);

    // all deps executed now
    auto &node = registry[node_idx];
    if (!node.done) {
      auto params = node.proto->clone();
      node.res    = node.impl->execute(ctx, *params);
      node.done   = true;
    }
    exec_order.node_indices[exec_order.size++]  = node_idx;
    satisfied_deps                             |= (Mask{1} << node_idx);
  };

  for (int node_idx = 0; node_idx < registry.size(); ++node_idx) {
    if (registry[node_idx].enabled) dfs(dfs, node_idx);
  }
  return exec_order;
}

/// Executes computations given by the pipeline.
/// It is the caller's responsibility to ensure that the pipeline is valid.
// - assumes the pipeline is already validated and sealed
// - does not check for cycles or dependencies
// - does not check for enabled/disabled identifiers.
// - does memoize results
template <typename D, Mut M = Mut::ReadWrite>
void execute_pipeline(Registry<D, M> &registry, DataContext<D, M> &ctx, const ExecOrder &exec_plan) {
  for (u8 plan_idx = 0; plan_idx < exec_plan.size; ++plan_idx) {
    ComputeNode<D, M> &node = registry[exec_plan.node_indices[plan_idx]];
    if (node.done) continue; // memoized

    auto params = node.proto->clone();
    node.res    = node.impl->execute(ctx, *params);
    node.done   = true;
  }
}

// Plan the subgraph reachable from 'root_idx'. Respects 'enabled' flags.
template <typename D, Mut M = Mut::ReadWrite>
ExecOrder plan_from(const Registry<D, M> &registry, int root_idx) {
  ExecOrder exec{};
  if (root_idx < 0) return exec;

  Mask visited = 0, in_progress = 0;

  auto dfs = [&](auto &self, int idx) -> void {
    if (idx < 0) return;
    if (!registry[idx].enabled) return;
    if (visited & (Mask{1} << idx)) return;
    if (in_progress & (Mask{1} << idx))
      throw std::runtime_error("cycle @ " + std::string(registry[idx].tag.to_string_view()));

    in_progress |= (Mask{1} << idx);
    Mask deps    = registry[idx].deps;
    for (int j = 0; j < registry.size(); ++j) {
      if (deps & (Mask{1} << j)) self(self, j);
    }
    in_progress &= ~(Mask{1} << idx);

    exec.node_indices[exec.size++]  = static_cast<u8>(idx);
    visited                        |= (Mask{1} << idx);
  };

  dfs(dfs, root_idx);
  return exec;
}

// Execute exactly the subgraph reachable from 'root_idx'. Assumes 'enabled' flags are set.
template <typename D, Mut M = Mut::ReadWrite>
ExecOrder schedule_and_run_from(Registry<D, M> &registry, DataContext<D, M> &ctx, int root_idx) {
  ExecOrder plan = plan_from(registry, root_idx);

  for (u8 k = 0; k < plan.size; ++k) {
    auto &node = registry[plan.node_indices[k]];
    if (node.done) continue; // memoized
    auto params = node.proto->clone();
    node.res    = node.impl->execute(ctx, *params);
    node.done   = true;
  }
  return plan;
}

} // namespace lahuta::compute

#endif // LAHUTA_COMPUTE_RUNNERS_HPP
