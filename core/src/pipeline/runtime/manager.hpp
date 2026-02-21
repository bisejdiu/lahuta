/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s = std::string{"besian"} + "sejdiu" + "@gmail.com";
 *   return *std::launder(&s);
 * }();
 *
 */

#ifndef LAHUTA_PIPELINE_RUNTIME_MANAGER_HPP
#define LAHUTA_PIPELINE_RUNTIME_MANAGER_HPP

#include <algorithm>
#include <atomic>
#include <chrono>
#include <functional>
#include <initializer_list>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "analysis/system/computation.hpp"
#include "analysis/topology/computation.hpp"
#include "logging/logging.hpp"
#include "pipeline/data/data_requirements.hpp"
#include "pipeline/ingest/adapters/lmdb.hpp"
#include "pipeline/ingest/descriptor.hpp"
#include "pipeline/ingest/realizer.hpp"
#include "pipeline/io/channel_multiplexer.hpp"
#include "pipeline/io/sink_iface.hpp"
#include "pipeline/metrics/run_metrics.hpp"
#include "pipeline/metrics/run_observer.hpp"
#include "pipeline/runtime/executor.hpp"
#include "pipeline/runtime/ingest_stream.hpp"
#include "pipeline/task/compute/computations.hpp"
#include "pipeline/task/compute/context.hpp"
#include "pipeline/task/compute/parameters.hpp"
#include "pipeline/task/task.hpp"
#include "runtime.hpp"

namespace lahuta::pipeline {
namespace C = lahuta::compute;
namespace P = lahuta::pipeline;

//
// StageManager orchestrates a DAG of named tasks over a stream of
// IngestDescriptors supplied by a source. Each descriptor (FileRef, LMDBRef,
// NMRRef, MDRef, ...) is turned into one or more PipelineItems by the
// DefaultRealizer. Items carry session metadata, conformer/frame ids, and an
// optional StreamSession/FrameHandle used to reuse heavyweight state across
// frames. For every emitted PipelineItem, the manager executes tasks in
// topological order with a fresh TaskContext and forwards Emissions to the
// configured sinks.
//
// Key concepts:
// - Descriptor sources: produce IngestDescriptor records describing where the
//   data lives. Helpers in ingest/ (Directory, Vector,
//   FileList, LMDBDescriptor, ...) expose next() so the
//   StageManager and other callers can feed the realizer uniformly.
// - DefaultRealizer: expands descriptors into PipelineItems. File/LMDB entries
//   emit a single item. Streaming descriptors (NMR/MD) emit one item per frame
//   while preserving StreamSession lifetime and backpressure limits.
// - Tasks (ITask): type-erased computations implementing run(path, ctx)
//   -> TaskResult. Dependencies declare the DAG that compile() validates.
// - Builtins: implicit "system" and "topology" computations injected when
//   requested, honoring parameter overrides.
// - TaskContext: per-item scratch space for typed objects and small strings.
//   Objects are shared as immutable shared_ptr instances between tasks.
// - Emissions & channels: TaskResult payloads are routed through the
//   ChannelMultiplexer so multiple tasks and sinks can share logical channels.
//
// Multi-run behavior & parameter invalidation:
// - run() may be called repeatedly. Sources are reset and the DefaultRealizer
//   is cleared before each run so descriptors can be reprocessed.
// - Parameter changes (via get_system_params/get_topology_params) trigger
//   invalidate_compilation(), forcing the compute registry to rebuild.
// - Each worker thread uses a ComputeEngine that is reset per item to avoid
//   leaking intermediate state across runs.
//
// Concurrency & safety:
// - Tasks register a thread_safe hint. If any task is marked unsafe the run is
//   forced to single-threaded execution. Otherwise the requested worker count
//   is used.
// - StreamSession::Permit acquired inside StageExecutor bounds the number of
//   in-flight frames per session, coordinating with data sources that require
//   serialized frame decoding.
// - TaskContext instances are confined to a worker while processing one item.
//   Emissions can interleave across items when multiple threads are active.
//
// Error handling:
// - A task returning ok=false short-circuits downstream tasks for that item but
//   leaves already emitted payloads untouched.
// - Dependency validation happens during compile(). Missing nodes or cycles
//   raise std::runtime_error.
//

class StageManager {
public:
  using SourcePtr = std::shared_ptr<IDescriptor>;

  enum class ReportingLevel { Off, Basic, Debug };

  struct StageTiming {
    std::string label;
    double setup_seconds   = 0.0;
    double compute_seconds = 0.0;
  };

  struct RunReport {
    // run identity & configuration
    std::size_t run_token         = 0;
    std::size_t stage_count       = 0;
    std::size_t threads_requested = 0;
    std::size_t threads_used      = 0;
    bool all_thread_safe          = true;
    bool metrics_enabled          = true;

    // item processing
    std::size_t items_total     = 0;
    std::size_t items_processed = 0;
    std::size_t items_skipped   = 0;

    // high-level timing
    double total_seconds = 0.0;
    double cpu_seconds   = 0.0;
    double io_seconds    = 0.0;

    // timing breakdown
    double ingest_seconds  = 0.0;
    double prepare_seconds = 0.0;
    double flush_seconds   = 0.0;
    double setup_seconds   = 0.0;
    double compute_seconds = 0.0;

    // per-stage breakdown (debuggin)
    std::vector<StageTiming> stage_breakdown;

    // queue/backpressure metrics
    std::size_t peak_inflight_items = 0;
    double average_queue_depth      = 0.0;

    // permit wait metrics
    std::size_t permit_wait_events   = 0;
    double permit_wait_total_seconds = 0.0;
    double permit_wait_min_seconds   = 0.0;
    double permit_wait_max_seconds   = 0.0;
    double permit_wait_avg_seconds   = 0.0;

    // multiplexer/sink I/O stats
    std::size_t mux_sink_count           = 0;
    std::uint64_t mux_enqueued_msgs      = 0;
    std::uint64_t mux_enqueued_bytes     = 0;
    std::uint64_t mux_written_msgs       = 0;
    std::uint64_t mux_written_bytes      = 0;
    std::uint64_t mux_stall_ns           = 0;
    std::uint64_t mux_drops              = 0;
    std::size_t mux_queue_depth_peak     = 0;
    std::size_t mux_queue_bytes_peak     = 0;
    std::size_t mux_active_writers_total = 0;
    std::size_t mux_active_writers_peak  = 0;
  };

  explicit StageManager(SourcePtr src) : src_(std::move(src)) {
    if (!src_) {
      throw std::invalid_argument("StageManager requires a descriptor source");
    }
    sys_params_.is_model = false;
    top_params_.flags    = TopologyComputation::All;
    Logger::get_logger()->debug("StageManager: created (auto_builtins={}, source_ptr={})",
                                auto_builtins_,
                                static_cast<const void *>(src_.get()));
  }

  // Builtin injection policy. Default false.
  void set_auto_builtins(bool on) {
    auto_builtins_ = on;
    invalidate_compilation();
    Logger::get_logger()->debug("StageManager: auto_builtins set to {}", auto_builtins_);
  }
  bool get_auto_builtins() const noexcept { return auto_builtins_; }

  void set_reporting_level(ReportingLevel level) { reporting_level_ = level; }
  ReportingLevel get_reporting_level() const noexcept { return reporting_level_; }

  void set_run_observer(std::shared_ptr<IRunObserver> observer) { run_observer_ = std::move(observer); }
  std::shared_ptr<IRunObserver> get_run_observer() const { return run_observer_; }

  void set_task_data_requirements(const std::string &name, P::DataFieldSet fields) {
    auto it = nodes_.find(name);
    if (it == nodes_.end()) {
      throw std::invalid_argument("StageManager: unknown task '" + name + "' in set_task_data_requirements");
    }
    it->second.requirements |= fields;
    invalidate_compilation();
  }

  P::DataFieldSet get_task_data_requirements(const std::string &name) const {
    auto it = nodes_.find(name);
    if (it == nodes_.end()) {
      throw std::invalid_argument("StageManager: unknown task '" + name + "' in get_task_data_requirements");
    }
    return it->second.requirements;
  }

  // Register or replace a task node with dependencies and implementation.
  // Internals:
  // - Records insertion order the first time a label is seen (in targets_)
  //   to drive deterministic per-item streaming order.
  // - Stores/updates the node in nodes_ and invalidates compute_factories_
  //   so compile() will rebuild the registry factories.
  // - thread_safe is used to decide if run() can parallelize items.
  void add_task(std::string name, std::vector<std::string> deps, std::shared_ptr<ITask> impl,
                bool thread_safe = true) {
    if (!impl) throw std::invalid_argument("StageManager.add_task: impl is null");
    bool fresh        = nodes_.find(name) == nodes_.end();
    std::string label = name;
    Node n;
    n.name        = std::move(name);
    n.deps        = std::move(deps);
    n.task        = std::move(impl);
    n.thread_safe = thread_safe;
    if (n.task) n.requirements = n.task->data_requirements();
    if (fresh) targets_.push_back(n.name);
    nodes_[n.name] = std::move(n);
    compute_factories_.clear(); // invalidate compiled factories
    Logger::get_logger()->debug("StageManager: added task '{}' (thread_safe={} fresh={})",
                                label,
                                thread_safe,
                                fresh);
  }

  // Subscribe a sink to a channel. Multiple sinks may subscribe to the same
  // channel. Multiple tasks may emit to the same channel.
  void connect_sink(const std::string &channel, std::shared_ptr<IDynamicSink> sink,
                    std::optional<BackpressureConfig> cfg = std::nullopt) {
    mux_.connect(channel, std::move(sink), std::move(cfg));
  }

  // Add a compute-backed task via a factory producing a C::Computation.
  // Internals:
  // - Records insertion order on first add (targets_) for streaming order.
  // - Stores/updates the node in nodes_ and invalidates compute_factories_.
  // Rationale:
  //   Using a factory decouples StageManager from specific computation types
  //   while letting compute introspect dependencies at compile-time.
  void add_computation(const std::string &name, std::vector<std::string> deps,
                       std::function<std::unique_ptr<C::Computation<P::PipelineContext>>()> factory,
                       bool thread_safe = true) {
    if (!factory) throw std::invalid_argument("StageManager.add_computation: factory is null");
    bool fresh        = nodes_.find(name) == nodes_.end();
    std::string label = name;
    Node n;
    n.name             = name;
    n.deps             = std::move(deps);
    n.task             = nullptr;
    n.thread_safe      = thread_safe;
    auto probe_factory = factory;
    try {
      auto sample = probe_factory();
      if (sample) n.requirements = sample->data_requirements();
    } catch (...) {
      n.requirements = P::DataFieldSet::none();
      throw;
    }
    n.make_compute = std::move(factory);
    if (fresh) targets_.push_back(n.name);
    nodes_[n.name] = std::move(n);
    compute_factories_.clear(); // invalidate compiled factories
    Logger::get_logger()->debug("StageManager: added computation '{}' (thread_safe={} fresh={})",
                                label,
                                thread_safe,
                                fresh);
  }

  // Prepare compute registry factories with conditional builtin injection.
  void compile() {
    auto logger = Logger::get_logger();
    logger->debug("StageManager: compile start (nodes={} auto_builtins={})", nodes_.size(), auto_builtins_);
    // Collect declared deps to decide builtin injection.
    auto collect_declared_deps = [&](const Node &n) -> std::vector<std::string> {
      std::unordered_set<std::string> out;
      for (const auto &d : n.deps)
        out.insert(d);
      if (n.make_compute) {
        try {
          auto c = n.make_compute();
          for (const auto &lbl : c->get_dependencies()) {
            out.insert(std::string(lbl.to_string_view()));
          }
        } catch (...) {
          throw std::runtime_error("StageManager: failed to inspect dependencies for node '" + n.name + "'");
        }
      }
      return std::vector<std::string>(out.begin(), out.end());
    };

    // Build a simple name into deps map, then evaluate builtin injection
    std::unordered_map<std::string, std::vector<std::string>> deps_by_name;
    for (auto &[name, node] : nodes_) {
      deps_by_name.emplace(name, collect_declared_deps(node));
    }

    const bool user_overrides_system   = nodes_.count("system") > 0;
    const bool user_overrides_topology = nodes_.count("topology") > 0;

    auto inj = decide_builtin_injection(auto_builtins_,
                                        deps_by_name,
                                        user_overrides_system,
                                        user_overrides_topology);

    // Persist factories
    compute_factories_.clear();
    compiled_builtins_.clear();

    // Injected builtins first (fixed labels and cheap construction)
    if (inj.inject_system) {
      compute_factories_.push_back(
          [this]() { return std::make_unique<analysis::SystemReadComputation>(sys_params_); });
      compiled_builtins_.insert("system");
      logger->debug("StageManager: injecting builtin 'system'");
    }
    if (inj.inject_topology) {
      compute_factories_.push_back(
          [this]() { return std::make_unique<analysis::BuildTopologyComputation>(top_params_); });
      compiled_builtins_.insert("topology");
      logger->debug("StageManager: injecting builtin 'topology'");
    }

    // User nodes (compute-backed or dynamic ITask). Order here is irrelevant.
    // Streaming uses targets_.
    for (const auto &[name, n] : nodes_) {
      if (n.make_compute) {
        compute_factories_.push_back(n.make_compute);
      } else if (n.task) {
        std::vector<std::string> sdeps = n.deps;
        auto task                      = n.task;
        compute_factories_.push_back([label = n.name, sdeps, task]() {
          return std::make_unique<DynamicTaskComputation>(label, sdeps, task);
        });
      }
    }

    graph_requirements_ = P::DataFieldSet::none();
    if (inj.inject_system) {
      auto sample          = std::make_unique<analysis::SystemReadComputation>(sys_params_);
      graph_requirements_ |= sample->data_requirements();
    }
    if (inj.inject_topology) {
      auto sample          = std::make_unique<analysis::BuildTopologyComputation>(top_params_);
      graph_requirements_ |= sample->data_requirements();
    }
    for (const auto &[name, node] : nodes_) {
      graph_requirements_ |= node.requirements;
    }
    realizer_.set_requirements(graph_requirements_);
    logger->debug("StageManager: compile finished (factories={} builtins={})",
                  compute_factories_.size(),
                  compiled_builtins_.size());
  }

  //
  // Execute the DAG over the configured source.
  // For each item: create a fresh TaskContext, run tasks in topo order, and
  // forward any Emissions to the ChannelMultiplexer. If any user-registered
  // node (compute-backed or dynamic ITask) is marked non-thread-safe, the run
  // executes single-threaded. Injected builtins are not considered in this
  // gate.
  //
  RunReport run(std::size_t threads = 1) {

    last_report_.reset();
    const auto total_begin = std::chrono::steady_clock::now();

    const std::size_t requested_threads = std::max<std::size_t>(threads, std::size_t{1});
    warn_if_lmdb_readers_exhausted(requested_threads);

    if (compute_factories_.empty()) compile();
    realizer_.set_requirements(graph_requirements_);

    LahutaRuntime::ensure_initialized(requested_threads); // Size thread dependent resources for the worker
                                                          // count

    // Source is reset so it is reusable across multiple run() calls.
    src_->reset();
    realizer_.reset();

    bool all_thread_safe = true;
    for (const auto &kv : nodes_) {
      all_thread_safe = all_thread_safe && kv.second.thread_safe;
    }

    //
    // Snapshot compiled plan and execute via StageExecutor
    // Snapshot validity: pointers reference manager-owned containers.
    // Between compile() and the end of this run(), these containers are stable.
    // Any parameter change triggers invalidate_compilation() before the next
    // run(), clearing/rebuilding containers so no stale reads will occur.
    //
    CompiledStage snapshot;
    snapshot.targets         = &targets_;
    snapshot.factories       = &compute_factories_;
    snapshot.all_thread_safe = all_thread_safe;
    snapshot.labels.clear();
    snapshot.labels.reserve(targets_.size());
    for (const auto &name : targets_) {
      snapshot.labels.emplace_back(name.c_str());
    }

    mux_.reopen_if_closed();
    const std::size_t run_token = 1 + global_run_epoch_.fetch_add(1, std::memory_order_relaxed);
    Logger::get_logger()->debug("StageManager: starting run (token={} "
                                "threads={} stages={} all_thread_safe={})",
                                run_token,
                                requested_threads,
                                snapshot.labels.size(),
                                all_thread_safe);

    StageMetricsSnapshot metrics_snapshot{};
    const bool metrics_enabled = reporting_level_ != ReportingLevel::Off;

    if (metrics_enabled) {
      StageRunMetrics metrics(reporting_level_ == ReportingLevel::Debug);
      metrics.configure_stage_breakdown(snapshot.labels.size());
      StageExecutor<StageRunMetrics> executor(snapshot,
                                              mux_,
                                              run_token,
                                              metrics,
                                              graph_requirements_,
                                              run_observer_);
      IngestItemStream item_stream(*src_, realizer_);
      executor.run(item_stream, requested_threads);

      const auto flush_begin = std::chrono::steady_clock::now();
      mux_.close_and_flush(flush_timeout_);
      const auto flush_end = std::chrono::steady_clock::now();
      metrics.add_flush(to_ns(flush_end - flush_begin));

      metrics_snapshot = metrics.snapshot();
    } else {
      NullStageRunMetrics metrics;
      StageExecutor<NullStageRunMetrics> executor(snapshot,
                                                  mux_,
                                                  run_token,
                                                  metrics,
                                                  graph_requirements_,
                                                  run_observer_);
      IngestItemStream item_stream(*src_, realizer_);
      executor.run(item_stream, requested_threads);

      const auto flush_begin = std::chrono::steady_clock::now();
      mux_.close_and_flush(flush_timeout_);
      const auto flush_end = std::chrono::steady_clock::now();
      metrics.add_flush(to_ns(flush_end - flush_begin));

      metrics_snapshot = metrics.snapshot();
    }

    const auto sink_stats = mux_.stats();

    const auto total_end = std::chrono::steady_clock::now();

    auto to_seconds = [](std::chrono::nanoseconds ns) { return std::chrono::duration<double>(ns).count(); };

    const auto total_ns = to_ns(total_end - total_begin);

    RunReport report;
    report.total_seconds   = to_seconds(total_ns);
    report.ingest_seconds  = to_seconds(std::chrono::nanoseconds(metrics_snapshot.ingest_ns));
    report.prepare_seconds = to_seconds(std::chrono::nanoseconds(metrics_snapshot.prepare_ns));
    report.flush_seconds   = to_seconds(std::chrono::nanoseconds(metrics_snapshot.flush_ns));
    report.setup_seconds   = to_seconds(std::chrono::nanoseconds(metrics_snapshot.setup_ns));
    report.compute_seconds = to_seconds(std::chrono::nanoseconds(metrics_snapshot.compute_ns));

    const auto cpu_ns  = metrics_snapshot.setup_ns + metrics_snapshot.compute_ns;
    const auto io_ns   = metrics_snapshot.ingest_ns + metrics_snapshot.prepare_ns + metrics_snapshot.flush_ns;
    report.cpu_seconds = to_seconds(std::chrono::nanoseconds(cpu_ns));
    report.io_seconds  = to_seconds(std::chrono::nanoseconds(io_ns));

    report.items_total         = metrics_snapshot.items_total;
    report.items_skipped       = metrics_snapshot.items_skipped;
    report.items_processed     = report.items_total >= report.items_skipped
                                     ? (report.items_total - report.items_skipped)
                                     : 0;
    report.stage_count         = snapshot.labels.size();
    report.threads_requested   = requested_threads;
    report.all_thread_safe     = all_thread_safe;
    report.threads_used        = all_thread_safe ? requested_threads : std::size_t{1};
    report.run_token           = run_token;
    report.metrics_enabled     = metrics_enabled;
    report.peak_inflight_items = metrics_snapshot.inflight_peak;
    if (metrics_snapshot.inflight_samples > 0) {
      report.average_queue_depth = static_cast<double>(metrics_snapshot.inflight_sum) /
                                   static_cast<double>(metrics_snapshot.inflight_samples);
    }
    report.permit_wait_events = static_cast<std::size_t>(metrics_snapshot.permit_wait_samples);

    // clang-format off
    if (report.permit_wait_events > 0) {
      report.permit_wait_total_seconds = to_seconds(std::chrono::nanoseconds(metrics_snapshot.permit_wait_ns_total));
      report.permit_wait_min_seconds   = to_seconds(std::chrono::nanoseconds(metrics_snapshot.permit_wait_ns_min));
      report.permit_wait_max_seconds   = to_seconds(std::chrono::nanoseconds(metrics_snapshot.permit_wait_ns_max));
      report.permit_wait_avg_seconds   = report.permit_wait_total_seconds /
                                         static_cast<double>(report.permit_wait_events);
    }
    // clang-format on
    if (metrics_enabled && reporting_level_ == ReportingLevel::Debug &&
        !metrics_snapshot.stage_compute_ns.empty()) {
      const auto stage_count = std::min({metrics_snapshot.stage_compute_ns.size(),
                                         metrics_snapshot.stage_setup_ns.size(),
                                         snapshot.labels.size()});
      report.stage_breakdown.reserve(stage_count);
      for (std::size_t i = 0; i < stage_count; ++i) {
        StageTiming timing;
        timing.label           = std::string(snapshot.labels[i].to_string_view());
        timing.setup_seconds   = to_seconds(std::chrono::nanoseconds(metrics_snapshot.stage_setup_ns[i]));
        timing.compute_seconds = to_seconds(std::chrono::nanoseconds(metrics_snapshot.stage_compute_ns[i]));
        report.stage_breakdown.push_back(std::move(timing));
      }
    }
    std::uint64_t mux_enqueued_msgs      = 0;
    std::uint64_t mux_enqueued_bytes     = 0;
    std::uint64_t mux_written_msgs       = 0;
    std::uint64_t mux_written_bytes      = 0;
    std::uint64_t mux_stall_ns           = 0;
    std::uint64_t mux_drops              = 0;
    std::size_t mux_queue_depth_peak     = 0;
    std::size_t mux_queue_bytes_peak     = 0;
    std::size_t mux_active_writers_total = 0;
    std::size_t mux_active_writers_peak  = 0;
    for (const auto &s : sink_stats) {
      mux_enqueued_msgs  += s.enqueued_msgs;
      mux_enqueued_bytes += s.enqueued_bytes;
      mux_written_msgs   += s.written_msgs;
      mux_written_bytes  += s.written_bytes;
      mux_stall_ns       += s.stalled_ns;
      mux_drops          += s.drops;
      if (s.queue_high_water_msgs > mux_queue_depth_peak) mux_queue_depth_peak = s.queue_high_water_msgs;
      if (s.queue_high_water_bytes > mux_queue_bytes_peak) mux_queue_bytes_peak = s.queue_high_water_bytes;
      mux_active_writers_total += s.writer_threads;
      if (s.writer_threads > mux_active_writers_peak) mux_active_writers_peak = s.writer_threads;
    }
    report.mux_sink_count           = sink_stats.size();
    report.mux_enqueued_msgs        = mux_enqueued_msgs;
    report.mux_enqueued_bytes       = mux_enqueued_bytes;
    report.mux_written_msgs         = mux_written_msgs;
    report.mux_written_bytes        = mux_written_bytes;
    report.mux_stall_ns             = mux_stall_ns;
    report.mux_drops                = mux_drops;
    report.mux_queue_depth_peak     = mux_queue_depth_peak;
    report.mux_queue_bytes_peak     = mux_queue_bytes_peak;
    report.mux_active_writers_total = mux_active_writers_total;
    report.mux_active_writers_peak  = mux_active_writers_peak;

    Logger::get_logger()->debug("StageManager: completed run_token {} (total={:.6f}s cpu={:.6f}s "
                                "io={:.6f}s items={} skipped={})",
                                run_token,
                                report.total_seconds,
                                report.cpu_seconds,
                                report.io_seconds,
                                report.items_total,
                                report.items_skipped);

    last_report_ = report;
    return report;
  }

  const std::vector<std::string> &sorted_tasks() const { return targets_; }
  const std::optional<RunReport> &last_report() const noexcept { return last_report_; }
  P::DataFieldSet data_requirements() const noexcept { return graph_requirements_; }

  // Expose current sink stats snapshot for diagnostics/observability.
  std::vector<ChannelMultiplexer::SinkStatsSnapshot> stats() const { return mux_.stats(); }

  // Flush timeout configuration
  void set_flush_timeout(std::chrono::milliseconds timeout) {
    if (timeout.count() < 0) {
      throw std::invalid_argument("StageManager.set_flush_timeout: timeout must be non-negative");
    }
    flush_timeout_ = timeout;
  }
  std::chrono::milliseconds get_flush_timeout() const { return flush_timeout_; }

  // Parameter access for builtins
  SystemReadParams &get_system_params() { return sys_params_; }
  const SystemReadParams &get_system_params() const { return sys_params_; }

  BuildTopologyParams &get_topology_params() { return top_params_; }
  const BuildTopologyParams &get_topology_params() const { return top_params_; }

  // Invalidate compiled state when parameters change
  void invalidate_compilation() {
    /*sorted_.clear();*/
    compute_factories_.clear();
    compiled_builtins_.clear();
    graph_requirements_ = P::DataFieldSet::none();
  }

  // Graph introspection for Python API
  std::vector<std::pair<std::string, std::vector<std::string>>> describe_graph() const {
    std::vector<std::pair<std::string, std::vector<std::string>>> result;
    // Only include builtins if present (injected or user-defined overrides)
    auto include_builtin = [&](const char *name, std::vector<std::string> deps) {
      if (compiled_builtins_.count(name) || nodes_.count(name)) {
        result.emplace_back(name, std::move(deps));
      }
    };
    include_builtin("system", {});
    include_builtin("topology", {"system"});
    for (const auto &[name, node] : nodes_) {
      result.emplace_back(name, node.deps);
    }
    return result;
  }

private:
  void warn_if_lmdb_readers_exhausted(std::size_t requested_threads) const {
    auto *lmdb_src = dynamic_cast<LMDBAdapter *>(src_.get());
    if (!lmdb_src) return;

    auto max_readers = lmdb_src->max_readers();
    if (!max_readers || *max_readers == 0) return;

    const std::size_t required_readers = requested_threads + 1; // here +1 for the main thread reading keys.
    if (required_readers > *max_readers) {
      throw std::runtime_error("StageManager: LMDB max_readers (" + std::to_string(*max_readers) +
                               ") is below required_readers (" + std::to_string(required_readers) +
                               ", threads + 1). Increase LAHUTA_LMDB_MAXREADERS, pass "
                               "LMDBEnvOptions "
                               "when constructing the LMDB source, or reduce the thread count.");
    }
  }

  struct BuiltinInjectionDecision {
    bool need_system;
    bool need_topology;
    bool inject_system;
    bool inject_topology;
  };

  static inline BuiltinInjectionDecision
  decide_builtin_injection(bool auto_builtins,
                           const std::unordered_map<std::string, std::vector<std::string>> &deps_by_name,
                           bool user_overrides_system, bool user_overrides_topology) {
    bool need_system   = false;
    bool need_topology = false;
    for (const auto &kv : deps_by_name) {
      for (const auto &d : kv.second) {
        if (d == "topology") need_topology = true;
        if (d == "system") need_system = true;
      }
    }
    if (need_topology) need_system = true; // topology implies system

    bool inject_system   = auto_builtins && need_system && !user_overrides_system;
    bool inject_topology = auto_builtins && need_topology && !user_overrides_topology;

    return BuiltinInjectionDecision{need_system, need_topology, inject_system, inject_topology};
  }

  struct Node {
    std::string name;
    std::vector<std::string> deps;
    std::shared_ptr<ITask> task;
    bool thread_safe = true;
    std::function<std::unique_ptr<C::Computation<P::PipelineContext>>()> make_compute;
    P::DataFieldSet requirements = P::DataFieldSet::none();
  };

private:
  SourcePtr src_;
  Realizer realizer_;
  std::unordered_map<std::string, Node> nodes_;
  std::vector<std::string> targets_;
  std::vector<std::function<std::unique_ptr<C::Computation<P::PipelineContext>>()>> compute_factories_;
  ChannelMultiplexer mux_;

  // Global run epoch: increments every run() call globally to scope TLS reuse
  // to a single run
  static inline std::atomic<std::size_t> global_run_epoch_{0};

  // builtin parameter storage
  SystemReadParams sys_params_;
  BuildTopologyParams top_params_;

  // Policy + introspection state
  bool auto_builtins_ = false;
  std::unordered_set<std::string> compiled_builtins_;
  ReportingLevel reporting_level_{ReportingLevel::Basic};
  std::shared_ptr<IRunObserver> run_observer_;
  P::DataFieldSet graph_requirements_ = P::DataFieldSet::none();

  // Flush timeout configuration
  std::chrono::milliseconds flush_timeout_{std::chrono::seconds(60)};

  std::optional<RunReport> last_report_;
};

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_RUNTIME_MANAGER_HPP
