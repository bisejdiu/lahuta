#ifndef LAHUTA_PIPELINE_DYNAMIC_MANAGER_HPP
#define LAHUTA_PIPELINE_DYNAMIC_MANAGER_HPP

#include <algorithm>
#include <atomic>
#include <chrono>
#include <functional>
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
#include "logging.hpp"
#include "pipeline/compute/computations.hpp"
#include "pipeline/compute/context.hpp"
#include "pipeline/compute/parameters.hpp"
#include "pipeline/dynamic/channel_multiplexer.hpp"
#include "pipeline/dynamic/executor.hpp"
#include "pipeline/dynamic/ingest_stream.hpp"
#include "pipeline/dynamic/sink_iface.hpp"
#include "pipeline/dynamic/types.hpp"
#include "runtime.hpp"
#include "sources/descriptor.hpp"
#include "sources/realizer.hpp"

// clang-format off
namespace lahuta::pipeline::dynamic {

template <class... Ts>
struct Overloaded : Ts... { using Ts::operator()...; };
template <class... Ts>
Overloaded(Ts...)->Overloaded<Ts...>;

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
//   data lives. Helpers in sources/ (Directory, Vector,
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
  using SourcePtr = std::shared_ptr<sources::IDescriptor>;

  struct RunReport {
    double total_seconds   = 0.0;
    double cpu_seconds     = 0.0;
    double io_seconds      = 0.0;
    double ingest_seconds  = 0.0;
    double prepare_seconds = 0.0;
    double flush_seconds   = 0.0;
    double setup_seconds   = 0.0;
    double compute_seconds = 0.0;

    std::size_t items_total     = 0;
    std::size_t items_processed = 0;
    std::size_t items_skipped   = 0;
    std::size_t stage_count     = 0;
    std::size_t threads_requested = 0;
    std::size_t threads_used      = 0;
    bool        all_thread_safe   = true;
    std::size_t run_token         = 0;
  };

  explicit StageManager(SourcePtr src)
      : src_(std::move(src)) {
    if (!src_) {
      throw std::invalid_argument("StageManager requires a descriptor source");
    }
    sys_params_.is_model = false;
    top_params_.flags = TopologyComputation::All;
    Logger::get_logger()->debug("StageManager: created (auto_builtins={}, source_ptr={})", auto_builtins_, static_cast<const void*>(src_.get()));
  }

  // Builtin injection policy. Default false.
  void set_auto_builtins(bool on) {
    auto_builtins_ = on;
    invalidate_compilation();
    Logger::get_logger()->debug("StageManager: auto_builtins set to {}", auto_builtins_);
  }
  bool get_auto_builtins() const noexcept { return auto_builtins_; }

  // Register or replace a task node with dependencies and implementation.
  // Internals:
  // - Records insertion order the first time a label is seen (in targets_)
  //   to drive deterministic per-item streaming order.
  // - Stores/updates the node in nodes_ and invalidates compute_factories_
  //   so compile() will rebuild the registry factories.
  // - thread_safe is used to decide if run() can parallelize items.
  void add_task(std::string name, std::vector<std::string> deps, std::shared_ptr<ITask> impl, bool thread_safe = true) {
    if (!impl) throw std::invalid_argument("StageManager.add_task: impl is null");
    bool fresh = nodes_.find(name) == nodes_.end();
    std::string label = name;
    Node n; n.name = std::move(name); n.deps = std::move(deps); n.task = std::move(impl); n.thread_safe = thread_safe;
    if (fresh) targets_.push_back(n.name);
    nodes_[n.name] = std::move(n);
    compute_factories_.clear(); // invalidate compiled factories
    Logger::get_logger()->debug("StageManager: added task '{}' (thread_safe={} fresh={})", label, thread_safe, fresh);
  }

  // Subscribe a sink to a channel. Multiple sinks may subscribe to the same
  // channel. Multiple tasks may emit to the same channel.
  void connect_sink(const std::string& channel,
                    std::shared_ptr<IDynamicSink> sink,
                    std::optional<BackpressureConfig> cfg = std::nullopt) {
    mux_.connect(channel, std::move(sink), std::move(cfg));
  }

  // Add a compute-backed task via a factory producing a compute::Computation.
  // Internals:
  // - Records insertion order on first add (targets_) for streaming order.
  // - Stores/updates the node in nodes_ and invalidates compute_factories_.
  // Rationale:
  //   Using a factory decouples StageManager from specific computation types
  //   while letting compute introspect dependencies at compile-time.
  void add_computation(const std::string& name,
                       std::vector<std::string> deps,
                       std::function<std::unique_ptr<compute::Computation<compute::PipelineContext, compute::Mut::ReadWrite>>()> factory,
                       bool thread_safe = true) {
    if (!factory) throw std::invalid_argument("StageManager.add_computation: factory is null");
    bool fresh = nodes_.find(name) == nodes_.end();
    std::string label = name;
    Node n; n.name = name; n.deps = std::move(deps); n.task = nullptr; n.thread_safe = thread_safe; n.make_compute = std::move(factory);
    if (fresh) targets_.push_back(n.name);
    nodes_[n.name] = std::move(n);
    compute_factories_.clear(); // invalidate compiled factories
    Logger::get_logger()->debug("StageManager: added computation '{}' (thread_safe={} fresh={})", label, thread_safe, fresh);
  }

  // Prepare compute registry factories with conditional builtin injection.
  void compile() {
    auto logger = Logger::get_logger();
    logger->debug("StageManager: compile start (nodes={} auto_builtins={})", nodes_.size(), auto_builtins_);
    // Collect declared deps to decide builtin injection.
    auto collect_declared_deps = [&](const Node& n) -> std::vector<std::string> {
      std::unordered_set<std::string> out;
      for (const auto& d : n.deps) out.insert(d);
      if (n.make_compute) {
        try {
          auto c = n.make_compute();
          for (const auto& lbl : c->get_dependencies()) {
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
    for (auto& [name, node] : nodes_) {
      deps_by_name.emplace(name, collect_declared_deps(node));
    }

    const bool user_overrides_system   = nodes_.count("system")   > 0;
    const bool user_overrides_topology = nodes_.count("topology") > 0;

    auto inj = decide_builtin_injection(auto_builtins_, deps_by_name, user_overrides_system, user_overrides_topology);

    // Persist factories
    compute_factories_.clear();
    compiled_builtins_.clear();

    // Injected builtins first (fixed labels and cheap construction)
    if (inj.inject_system) {
      compute_factories_.push_back([this]() {
        return std::make_unique<analysis::system::SystemReadComputation>(sys_params_);
      });
      compiled_builtins_.insert("system");
      logger->debug("StageManager: injecting builtin 'system'");
    }
    if (inj.inject_topology) {
      compute_factories_.push_back([this]() {
        return std::make_unique<analysis::topology::BuildTopologyComputation>(top_params_);
      });
      compiled_builtins_.insert("topology");
      logger->debug("StageManager: injecting builtin 'topology'");
    }

    // User nodes (compute-backed or dynamic ITask). Order here is irrelevant. Streaming uses targets_.
    for (const auto& [name, n] : nodes_) {
      if (n.make_compute) {
        compute_factories_.push_back(n.make_compute);
      } else if (n.task) {
        std::vector<std::string> sdeps = n.deps;
        auto task = n.task;
        compute_factories_.push_back([label = n.name, sdeps, task]() {
          return std::make_unique<pipeline::compute::DynamicTaskComputation>(label, sdeps, task);
        });
      }
    }
    logger->debug("StageManager: compile finished (factories={} builtins={})", compute_factories_.size(), compiled_builtins_.size());
  }

  //
  // Execute the DAG over the configured source.
  // For each item: create a fresh TaskContext, run tasks in topo order, and
  // forward any Emissions to the ChannelMultiplexer. If any user-registered node
  // (compute-backed or dynamic ITask) is marked non-thread-safe, the run executes
  // single-threaded. Injected builtins are not considered in this gate.
  //
  RunReport run(std::size_t threads = 1) {

    last_report_.reset();
    const auto total_begin = std::chrono::steady_clock::now();

    const std::size_t requested_threads = std::max<std::size_t>(threads, std::size_t{1});

    if (compute_factories_.empty()) compile();

    LahutaRuntime::ensure_initialized(requested_threads); // Size thread dependent resources for the worker count

    // Source is reset so it is reusable across multiple run() calls.
    src_->reset();
    realizer_.reset();

    bool all_thread_safe = true;
    for (const auto& kv : nodes_) all_thread_safe = all_thread_safe && kv.second.thread_safe;

    //
    // Snapshot compiled plan and execute via StageExecutor
    // Snapshot validity: pointers reference manager-owned containers.
    // Between compile() and the end of this run(), these containers are stable.
    // Any parameter change triggers invalidate_compilation() before the next run(),
    // clearing/rebuilding containers so no stale reads will occur.
    //
    CompiledStage snapshot;
    snapshot.targets   = &targets_;
    snapshot.factories = &compute_factories_;
    snapshot.all_thread_safe = all_thread_safe;
    snapshot.labels.clear();
    snapshot.labels.reserve(targets_.size());
    for (const auto& name : targets_) snapshot.labels.emplace_back(name.c_str());

    mux_.reopen_if_closed();
    const std::size_t run_token = 1 + global_run_epoch_.fetch_add(1, std::memory_order_relaxed);
    Logger::get_logger()->debug(
      "StageManager: starting run (token={} threads={} stages={} all_thread_safe={})",
      run_token, requested_threads, snapshot.labels.size(), all_thread_safe);

    StageExecutor::RunMetrics metrics;
    StageExecutor executor(snapshot, mux_, run_token, metrics);
    IngestItemStream item_stream(*src_, realizer_);
    executor.run(item_stream, requested_threads);

    // Stop ingress and drain writers with a deadline
    const auto flush_begin = std::chrono::steady_clock::now();
    mux_.close_and_flush(flush_timeout_);
    const auto flush_end = std::chrono::steady_clock::now();
    metrics.add_flush(std::chrono::duration_cast<std::chrono::nanoseconds>(flush_end - flush_begin));

    const auto total_end = std::chrono::steady_clock::now();

    auto to_seconds = [](std::chrono::nanoseconds ns) {
      return std::chrono::duration<double>(ns).count();
    };

    const auto total_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(total_end - total_begin);

    RunReport report;
    report.total_seconds     = to_seconds(total_ns);
    report.cpu_seconds       = to_seconds(metrics.cpu_duration());
    report.io_seconds        = to_seconds(metrics.io_duration());
    report.ingest_seconds    = to_seconds(metrics.ingest_duration());
    report.prepare_seconds   = to_seconds(metrics.prepare_duration());
    report.flush_seconds     = to_seconds(metrics.flush_duration());
    report.setup_seconds     = to_seconds(metrics.setup_duration());
    report.compute_seconds   = to_seconds(metrics.compute_duration());
    report.items_total       = metrics.items_total();
    report.items_skipped     = metrics.items_skipped();
    report.items_processed   = report.items_total >= report.items_skipped ? (report.items_total - report.items_skipped) : 0;
    report.stage_count       = snapshot.labels.size();
    report.threads_requested = requested_threads;
    report.all_thread_safe   = all_thread_safe;
    report.threads_used      = all_thread_safe ? requested_threads : std::size_t{1};
    report.run_token         = run_token;

    Logger::get_logger()->debug(
      "StageManager: completed run_token {} (total={:.6f}s cpu={:.6f}s io={:.6f}s items={} skipped={})",
      run_token,
      report.total_seconds,
      report.cpu_seconds,
      report.io_seconds,
      report.items_total,
      report.items_skipped);

    last_report_ = report;
    return report;
  }

  const std::vector<std::string>& sorted_tasks() const { return targets_; }

  const std::optional<RunReport>& last_report() const noexcept { return last_report_; }

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
  compute::SystemReadParams& get_system_params() { return sys_params_; }
  const compute::SystemReadParams& get_system_params() const { return sys_params_; }

  compute::BuildTopologyParams& get_topology_params() { return top_params_; }
  const compute::BuildTopologyParams& get_topology_params() const { return top_params_; }

  // Invalidate compiled state when parameters change
  void invalidate_compilation() {
    /*sorted_.clear();*/
    compute_factories_.clear();
    compiled_builtins_.clear();
  }

  // Graph introspection for Python API
  std::vector<std::pair<std::string, std::vector<std::string>>> describe_graph() const {
    std::vector<std::pair<std::string, std::vector<std::string>>> result;
    // Only include builtins if present (injected or user-defined overrides)
    auto include_builtin = [&](const char* name, std::vector<std::string> deps) {
      if (compiled_builtins_.count(name) || nodes_.count(name)) {
        result.emplace_back(name, std::move(deps));
      }
    };
    include_builtin("system",   {});
    include_builtin("topology", {"system"});
    for (const auto& [name, node] : nodes_) {
      result.emplace_back(name, node.deps);
    }
    return result;
  }

private:
  struct BuiltinInjectionDecision {
    bool need_system;
    bool need_topology;
    bool inject_system;
    bool inject_topology;
  };

  static inline BuiltinInjectionDecision decide_builtin_injection(
      bool auto_builtins,
      const std::unordered_map<std::string, std::vector<std::string>>& deps_by_name,
      bool user_overrides_system,
      bool user_overrides_topology) {
    bool need_system = false;
    bool need_topology = false;
    for (const auto& kv : deps_by_name) {
      for (const auto& d : kv.second) {
        if (d == "topology") need_topology = true;
        if (d == "system")   need_system   = true;
      }
    }
    if (need_topology) need_system = true; // topology implies system

    bool inject_system   = auto_builtins && need_system   && !user_overrides_system;
    bool inject_topology = auto_builtins && need_topology && !user_overrides_topology;

    return BuiltinInjectionDecision{need_system, need_topology, inject_system, inject_topology};
  }

  struct Node {
    std::string name;
    std::vector<std::string> deps;
    std::shared_ptr<ITask> task;
    bool thread_safe = true;
    std::function<std::unique_ptr<compute::Computation<compute::PipelineContext, compute::Mut::ReadWrite>>()> make_compute;
  };

private:
  SourcePtr src_;
  sources::Realizer realizer_;
  std::unordered_map<std::string, Node> nodes_;
  std::vector<std::string> targets_;
  std::vector<std::function<std::unique_ptr<compute::Computation<compute::PipelineContext, compute::Mut::ReadWrite>>()>> compute_factories_;
  ChannelMultiplexer mux_;

  // Global run epoch: increments every run() call globally to scope TLS reuse to a single run
  static inline std::atomic<std::size_t> global_run_epoch_{0};

  // builtin parameter storage
  pipeline::compute::SystemReadParams sys_params_;
  pipeline::compute::BuildTopologyParams top_params_;

  // Policy + introspection state
  bool auto_builtins_ = false;
  std::unordered_set<std::string> compiled_builtins_;

  // Flush timeout configuration
  std::chrono::milliseconds flush_timeout_{std::chrono::seconds(60)};

  std::optional<RunReport> last_report_;
};

} // namespace lahuta::pipeline::dynamic

#endif // LAHUTA_PIPELINE_DYNAMIC_MANAGER_HPP
