#ifndef LAHUTA_PIPELINE_DYNAMIC_MANAGER_HPP
#define LAHUTA_PIPELINE_DYNAMIC_MANAGER_HPP

#include <atomic>
#include <chrono>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <variant>
#include <vector>

#include "analysis/system/computation.hpp"
#include "analysis/topology/computation.hpp"
#include "db/db.hpp"
#include "runtime.hpp"
#include "pipeline/compute/computations.hpp"
#include "pipeline/compute/context.hpp"
#include "pipeline/compute/parameters.hpp"
#include "pipeline/core/emitter.hpp"
#include "pipeline/core/stage.hpp"
#include "pipeline/dynamic/types.hpp"
#include "pipeline/engine.hpp"
#include "sources/db_key_source.hpp"
#include "sources/directory_source.hpp"
#include "sources/file_list_source.hpp"
#include "sources/vector_source.hpp"
#include <pipeline/dynamic/channel_multiplexer.hpp>
#include <pipeline/dynamic/executor.hpp>
#include <pipeline/dynamic/sink_iface.hpp>

// clang-format off
namespace lahuta::pipeline::dynamic {

//
// StageManager orchestrates a DAG of named tasks over a configured Source.
// For each input item (e.g, a file path), tasks run in topological order with a fresh
// TaskContext. Tasks may write intermediate values into the context and/or
// return Emission records that are routed by channel to registered sinks.
//
// Key concepts:
// - Tasks (ITask): type-erased computations implementing run(path, ctx) -> TaskResult.
// - Dependencies (DAG): tasks declare upstream names. StageManager enforces a
//   topological execution order (compile() throws on missing deps or cycles).
// - builtins: implicit "system" and "topology" tasks with configurable parameters.
//   These are automatically included and managed by the StageManager.
// - TaskContext: per-item scratch space for typed objects and small strings.
//   Cross-task access is read-only for objects. Tasks update state by writing
//   new values into the context. Fresh context created for each item.
// - Emissions & channels: each TaskResult may contain Emission{channel,payload}.
//   The ChannelMultiplexer fans out per channel to all subscribed sinks
//   (many-to-many: multiple tasks -> same channel. Multiple sinks <- same channel).
//
// Multi-run behavior & parameter invalidation:
// - run() can be called multiple times on the same StageManager instance.
// - Parameter changes (via get_system_params/get_topology_params) trigger
//   invalidate_compilation(), forcing recompilation with updated parameters.
// - Sources are automatically reset before each run() to enable re-iteration
//   over the same input with different parameters.
// - Fresh ComputeEngine instances are created for each item in each run,
//   making sure no stale computation state persists across runs.
//
// Concurrency & safety:
// - User-registered nodes (compute-backed via add_computation or dynamic ITask via
//   add_task) carry a thread_safe hint. run() executes with a single worker if any
//   such node is marked unsafe, otherwise it uses the provided thread count.
// - Injected builtins ("system", "topology") are not considered in this gate and
//   are assumed thread-safe.
// - Python callables added through bindings preserve item-level parallelism: they
//   are registered as thread-safe for the stage. Optional per-task serialization is
//   handled inside the callable.
// - A given TaskContext is used by one worker for the duration of a single item.
// - Emissions are forwarded to the multiplexer as they are produced. Within an
//   item, ordering follows task order. Across items, emissions may interleave.
//
// Error handling:
// - If a task returns ok=false for an item, remaining downstream tasks for that
//   item are skipped. Emissions already produced for the item are still
//   delivered to sinks. No rollback is attempted.
// - Parameter validation occurs during compilation. Invalid dependencies or
//   cycles result in std::runtime_error.
//

template <typename T, typename = void>
struct has_reset : std::false_type {};

template <typename T>
struct has_reset<T, std::void_t<decltype(std::declval<T&>().reset())>> : std::true_type {};

template <typename T>
inline void reset_if_supported(T& obj) {
  if constexpr (has_reset<T>::value) {
    obj.reset();
  }
}

// Manages LMDBDatabase lifetime for a DB-backed source.
// Stores a shared_ptr<LMDBDatabase> and a DBKeySource referencing it.
struct DBKeySourceHolder {
  using value_type = std::string;

  // Construct by opening an LMDB environment at path with optional batch size
  explicit DBKeySourceHolder(const std::string& db_path, std::size_t batch_size = 1024)
      : db_(std::make_shared<lahuta::LMDBDatabase>(db_path)), src_(*db_, batch_size) {}

  // If a database handle is already available
  explicit DBKeySourceHolder(std::shared_ptr<lahuta::LMDBDatabase> db, std::size_t batch_size = 1024)
      : db_(std::move(db)), src_(*db_, batch_size) {}

  std::optional<std::string> next() { return src_.next(); }
  void reset() { src_.reset(); }

private:
  std::shared_ptr<lahuta::LMDBDatabase> db_;
  sources::DBKeySource src_;
};

class StageManager {
public:
  using SourceVariant = std::variant<
    sources::DirectorySource,
    sources::VectorSource,
    sources::FileListSource,
    DBKeySourceHolder
  >;

  explicit StageManager(SourceVariant src) : src_(std::move(src)) {
    // defaults
    sys_params_.is_model = false;
    top_params_.flags    = TopologyComputation::All;
  }

  // Host-agnostic policy knob. Python turns this on.
  void set_auto_builtins(bool on) {
    auto_builtins_ = on;
    invalidate_compilation();
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
    Node n; n.name = std::move(name); n.deps = std::move(deps); n.task = std::move(impl); n.thread_safe = thread_safe;
    if (fresh) targets_.push_back(n.name);
    nodes_[n.name] = std::move(n);
    compute_factories_.clear(); // invalidate compiled factories
  }

  // Subscribe a sink to a channel. Multiple sinks may subscribe to the same
  // channel. Multiple tasks may emit to the same channel.
  void connect_sink(const std::string& channel, std::shared_ptr<IDynamicSink> sink) {
    mux_.connect(channel, std::move(sink));
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
    Node n; n.name = name; n.deps = std::move(deps); n.task = nullptr; n.thread_safe = thread_safe; n.make_compute = std::move(factory);
    if (fresh) targets_.push_back(n.name);
    nodes_[n.name] = std::move(n);
    compute_factories_.clear(); // invalidate compiled factories
  }

  // Prepare compute registry factories with conditional builtin injection.
  void compile() {
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

    // Decide builtin injection by scanning declared deps
    struct GInfo { std::vector<std::string> deps; bool injected_builtin = false; };
    std::unordered_map<std::string, GInfo> G;
    for (auto& [name, node] : nodes_) {
      auto deps = collect_declared_deps(node);
      G.emplace(name, GInfo{std::move(deps), /*injected_builtin=*/false});
    }

    bool need_topology = false;
    bool need_system   = false;
    for (const auto& [name, gi] : G) {
      for (const auto& d : gi.deps) {
        if (d == "topology") need_topology = true;
        if (d == "system")   need_system   = true;
      }
    }
    if (need_topology) need_system = true; // topology -> system

    const bool user_overrides_system   = nodes_.count("system")   > 0;
    const bool user_overrides_topology = nodes_.count("topology") > 0;

    const bool inject_system   = auto_builtins_ && need_system   && !user_overrides_system;
    const bool inject_topology = auto_builtins_ && need_topology && !user_overrides_topology;

    // Persist factories
    compute_factories_.clear();
    compiled_builtins_.clear();

    // Injected builtins first (labels are fixed, construction is cheap)
    if (inject_system) {
      compute_factories_.push_back([this]() {
        return std::make_unique<analysis::system::SystemReadComputation>(sys_params_);
      });
      compiled_builtins_.insert("system");
    }
    if (inject_topology) {
      compute_factories_.push_back([this]() {
        return std::make_unique<analysis::topology::BuildTopologyComputation>(top_params_);
      });
      compiled_builtins_.insert("topology");
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
  }

  // Execute the DAG over the configured source.
  // For each item: create a fresh TaskContext, run tasks in topo order, and
  // forward any Emissions to the ChannelMultiplexer. If any user-registered node
  // (compute-backed or dynamic ITask) is marked non-thread-safe, the run executes
  // single-threaded. Injected builtins are not considered in this gate.
  void run(std::size_t threads = 1) {
    if (compute_factories_.empty()) {
      compile();
    }

    // Size thread dependent resources for the worker count
    lahuta::LahutaRuntime::ensure_initialized(threads);

    // Reset source to make it reusable across multiple run() calls.
    std::visit([](auto& src) { reset_if_supported(src); }, src_);

    bool all_thread_safe = true;
    for (const auto& kv : nodes_) all_thread_safe = all_thread_safe && kv.second.thread_safe;

    // Snapshot compiled plan and execute via StageExecutor
    // Snapshot validity: pointers reference manager-owned containers. Between
    // compile() and the end of this run(), these containers are stable. Any
    // parameter change triggers invalidate_compilation() before the next run(),
    // clearing/rebuilding containers so no stale reads will occur.
    CompiledStage snapshot;
    snapshot.targets   = &targets_;
    snapshot.factories = &compute_factories_;
    snapshot.all_thread_safe = all_thread_safe;
    snapshot.labels.clear();
    snapshot.labels.reserve(targets_.size());
    for (const auto& name : targets_) snapshot.labels.emplace_back(name.c_str());

    mux_.reopen_if_closed();
    const std::size_t run_token = 1 + global_run_epoch_.fetch_add(1, std::memory_order_relaxed);
    StageExecutor executor(snapshot, mux_, run_token);
    executor.run(src_, threads);

    // Stop ingress and drain writers with a deadline
    mux_.close_and_flush(std::chrono::seconds(5));
  }

  const std::vector<std::string>& sorted_tasks() const { return targets_; }

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
  struct Node {
    std::string name;
    std::vector<std::string> deps;
    std::shared_ptr<ITask> task;
    bool thread_safe = true;
    std::function<std::unique_ptr<compute::Computation<compute::PipelineContext, compute::Mut::ReadWrite>>()> make_compute;
  };

  template<typename Src>
  void run_on_source(Src& src, pipeline::Stage<std::string, void>& st, pipeline::IEmitter<void>& out, std::size_t threads) {
    pipeline::PipelineEngine{threads}.run(src, st, out);
  }

  void dsl_run(pipeline::Stage<std::string, void>& st, pipeline::IEmitter<void>& out, std::size_t threads) {
    std::visit([&](auto& src) {
      run_on_source(src, st, out, threads);
    }, src_);
  }

private:
  SourceVariant src_;
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
};

} // namespace lahuta::pipeline::dynamic

#endif // LAHUTA_PIPELINE_DYNAMIC_MANAGER_HPP
