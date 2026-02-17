/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s(22, 'X');
 *   auto it = s.begin();
 *   for (std::string_view part : {"besian", "sejdiu", "@gmail.com"}) {
 *     it = std::copy(part.begin(), part.end(), it);
 *   }
 *   return s;
 * }();
 *
 */

#ifndef LAHUTA_COMPUTE_ENGINE_HPP
#define LAHUTA_COMPUTE_ENGINE_HPP

#include <algorithm>
#include <cassert>
#include <deque>
#include <memory>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include "node.hpp"
#include "plan.hpp"
#include "registry.hpp"
#include "runners.hpp"

namespace lahuta::compute {
namespace C = lahuta::compute;

struct ResultEntry {
  ComputationLabel label;
  ComputationResult result;
};

template <typename DataT, Mut M>
class ComputeEngine { // default for M is set in context.hpp forward declaration
private:
  DataContext<DataT, M> ctx;
  Registry<DataT, M> registry;
  bool auto_heal_ = true;

  void run_impl(ComputationLabel root) {
    registry.seal();                 // idempotent, does nothing on 2nd call
    enable_chain(root);              // auto-enable full dependency tree
    schedule_and_run(registry, ctx); // executes + memoises + cycle-check
  }

  void enable_chain(ComputationLabel root) {

    if (!auto_heal_) return;

    std::deque<int> queue{registry.find(root)};
    Mask visited = 0;

    while (!queue.empty()) {
      int node_idx = queue.front();
      queue.pop_front();
      if (visited & (Mask{1} << node_idx)) continue;
      visited |= (Mask{1} << node_idx);

      if (!registry[node_idx].forced_disabled) {
        if (!registry[node_idx].enabled) { // resurrect disabled nodes
          registry[node_idx].enabled = true;
          invalidate_downstream(node_idx); // invalidate cached dependents
        }
      }

      Mask dep_mask = registry[node_idx].deps;
      for (int dep_idx = 0; dep_idx < registry.size(); ++dep_idx)
        if (dep_mask & (Mask{1} << dep_idx)) queue.push_back(dep_idx);
    }
  }

  void invalidate_downstream(int idx) {
    std::deque<int> queue{idx};
    Mask visited = 0;

    while (!queue.empty()) {
      int node_idx = queue.front();
      queue.pop_front();
      if (visited & (Mask{1} << node_idx)) continue;
      visited |= (Mask{1} << node_idx);

      registry[node_idx].done = false; // invalidate cached result

      Mask rev_dep_mask = registry[node_idx].rdeps; // nodes that depend on this node
      for (int dep_idx = 0; dep_idx < registry.size(); ++dep_idx)
        if (rev_dep_mask & (Mask{1} << dep_idx)) queue.push_back(dep_idx);
    }
  }

public:
  explicit ComputeEngine(DataT &data) : ctx(data, /*engine=*/this) {}

  /// Add a computation to the registry
  void add(std::unique_ptr<Computation<DataT, M>> c) { registry.add(std::move(c)); }

  /// If true, the engine will automatically enable dependent computations even if they are disabled.
  void set_auto_heal(bool on) { auto_heal_ = on; }

  /// value-returning computations  (R != void) -> returns std::optional<R>
  template <typename R, std::enable_if_t<!std::is_void_v<R>, int> = 0>
  std::optional<R> run(ComputationLabel root) { // target label
    run_impl(root);
    auto &node = registry[registry.find(root)];
    if (!node.done || !node.res.is_success()) return std::nullopt;
    return node.res.template get_value<R>();
  }

  /// "fire-and-forget" computations  (R == void) returns bool
  template <typename R = void, std::enable_if_t<std::is_void_v<R>, int> = 0>
  bool run(ComputationLabel root) { // target label
    run_impl(root);
    auto &node = registry[registry.find(root)];
    return node.done && node.res.is_success();
  }

  /// Execute a precomputed plan without re-planning.
  void execute_plan(const ExecOrder &plan) { //
    execute_pipeline(registry, ctx, plan);
  }

  void enable(ComputationLabel label, bool on) {
    int idx = registry.find(label);
    if (idx < 0) throw std::runtime_error("unknown label");

    registry[idx].enabled = on;
    /*reg_[i].forced_disabled = !on;*/
    invalidate_downstream(idx);
  }

  // Plan the subgraph from a root label (no execution). Applies auto-heal.
  ExecOrder plan_from(ComputationLabel root) {
    registry.seal();
    enable_chain(root);
    int idx = registry.find(root);
    if (idx < 0) throw std::runtime_error("unknown label");
    return C::plan_from(registry, idx);
  }

  // Execute exactly the subgraph reachable from 'root'. Applies auto-heal.
  template <typename R = void>
  bool run_from(ComputationLabel root) {
    registry.seal();
    enable_chain(root);
    int idx = registry.find(root);
    if (idx < 0) throw std::runtime_error("unknown label");
    (void)C::schedule_and_run_from(registry, ctx, idx);
    auto &node = registry[idx];
    return node.done && node.res.is_success();
  }

  /// Returns the result of a computation
  template <typename R>
  std::optional<R> result(ComputationLabel lbl) const {
    int idx = registry.find(lbl);
    if (idx < 0) return std::nullopt;
    auto &node = registry[idx];
    if (!node.done || !node.res.is_success()) return std::nullopt;

    try {
      return node.res.template get_value<R>();
    } catch (...) {
      return std::nullopt;
    }
  }

  /// Get the result of a computation by label. Triggers on-complete hook once.
  ComputationResult get_computation_result(ComputationLabel label) {
    int idx = registry.find(label);
    if (idx < 0 || !registry[idx].done) {
      return ComputationResult(ComputationError("Result not available"));
    }
    // TODO: Is the on_complete hook worth having?
    // Invoke per-computation on_complete exactly once
    auto &node = registry[idx];
    if (!node.postprocessed && node.impl) {
      try {
        node.impl->on_complete(ctx, node.res);
        // hook failing is not an issue
      } catch (...) {
      }
      node.postprocessed = true;
    }
    return node.res;
  }

  /// Checks if a computation has completed successfully
  bool has_completed(ComputationLabel label) const {
    int idx = registry.find(label);
    return (idx >= 0 && registry[idx].done && registry[idx].res.is_success());
  }

#ifdef LAHUTA_TESTING
  /// Test-only: execution count for a computation label
  std::uint64_t get_run_count(ComputationLabel label) const {
    int idx = registry.find(label);
    if (idx < 0) throw std::runtime_error("unknown label");
    return registry[idx].run_count;
  }

  /// Test-only: reset all execution counters
  void reset_run_counts() {
    for (int node_idx = 0; node_idx < registry.size(); ++node_idx) {
      registry[node_idx].run_count = 0;
    }
  }
#endif

  /// Returns all results of enabled computations
  std::vector<ResultEntry> w_execute_all() {
    registry.seal();

    // auto-heal the entire graph once
    if (auto_heal_) {
      for (int node_idx = 0; node_idx < registry.size(); ++node_idx)
        if (registry[node_idx].enabled) {
          enable_chain(registry[node_idx].tag); // may enable more nodes
        }
    }

    // run every enabled node
    for (int node_idx = 0; node_idx < registry.size(); ++node_idx)
      if (registry[node_idx].enabled) {
        run_impl(registry[node_idx].tag);
      }

    // move results out
    std::vector<ResultEntry> result;
    for (int node_idx = 0; node_idx < registry.size(); ++node_idx)
      if (registry[node_idx].enabled) {
        result.push_back({registry[node_idx].tag, registry[node_idx].res});
      }
    return result;
  }

  /// Execute all computations and move results out
  std::vector<ResultEntry> execute_and_move_all() {
    // seal & auto-heal exactly like execute_all()
    registry.seal();

    if (auto_heal_) {
      for (int node_idx = 0; node_idx < registry.size(); ++node_idx)
        if (registry[node_idx].enabled) {
          enable_chain(registry[node_idx].tag);
        }
    }

    // run every enabled node
    for (int node_idx = 0; node_idx < registry.size(); ++node_idx)
      if (registry[node_idx].enabled) {
        run_impl(registry[node_idx].tag);
      }

    // move results out + reset engine state
    std::vector<ResultEntry> result;
    result.reserve(registry.size());

    for (int node_idx = 0; node_idx < registry.size(); ++node_idx) {
      if (!registry[node_idx].enabled) continue;

      result.push_back({registry[node_idx].tag, std::move(registry[node_idx].res)});

      // clear cached state so next run recomputes
      registry[node_idx].res  = {};
      registry[node_idx].done = false;
    }
    return result;
  }

  /// drop cached results
  void reset() {
    for (int node_idx = 0; node_idx < registry.size(); ++node_idx) {
      registry[node_idx].res  = {};
      registry[node_idx].done = false;
    }
  }

  /// get parameters of a computation
  template <typename P>
  P &get_parameters(ComputationLabel label) {
    int node_idx = registry.find(label);
    if (node_idx < 0) throw std::runtime_error("unknown label");

    auto *base = registry[node_idx].proto.get();
    if (!base) throw std::runtime_error("no parameters stored");
    auto *typed = dynamic_cast<P *>(base);
    if (!typed) throw std::runtime_error("wrong parameter type");

    invalidate_downstream(node_idx); // invalidate cached dependents

    return *typed; // caller mutates in-place
  }

  /// run a pipeline of computations
  template <const ComputationLabel &...Order>
  void run_pipeline() {
    registry.seal();
    ExecutionPlan<DataT, M> exec_plan(registry, {Order...});
    execute_pipeline(registry, ctx, exec_plan.plan());
  }

  int find_label(const ComputationLabel &label) const noexcept { //
    return registry.find(label);
  }

  bool is_enabled(int idx) const noexcept {
    return idx >= 0 && idx < registry.size() && registry[idx].enabled;
  }

  /// check if a computation is available, i.e. enabled, returns true if yes
  bool is_computation_available(const ComputationLabel &lbl) const noexcept {
    int idx = find_label(lbl);
    return is_enabled(idx);
  }
};

} // namespace lahuta::compute

#endif // LAHUTA_COMPUTE_ENGINE_HPP
