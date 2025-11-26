#ifndef LAHUTA_PIPELINE_DYNAMIC_EXECUTOR_HPP
#define LAHUTA_PIPELINE_DYNAMIC_EXECUTOR_HPP

#include <atomic>
#include <chrono>
#include <cstddef>
#include <functional>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <GraphMol/Conformer.h>
#include <rdkit/Geometry/point.h>

#include "compute/engine.hpp"
#include "logging.hpp"
#include "pipeline/compute/context.hpp"
#include "pipeline/core/emitter.hpp"
#include "pipeline/core/stage.hpp"
#include "pipeline/data_requirements.hpp"
#include "pipeline/dynamic/channel_multiplexer.hpp"
#include "pipeline/dynamic/keys.hpp"
#include "pipeline/dynamic/run_observer.hpp"
#include "pipeline/dynamic/types.hpp"
#include "pipeline/engine.hpp"
#include "pipeline/model_payload.hpp"
#include "pipeline/pipeline_item.hpp"
#include "pipeline/stream_session.hpp"

// clang-format off
namespace lahuta::pipeline::dynamic {
using namespace lahuta::topology::compute;

template <class Dur>
constexpr std::chrono::nanoseconds to_ns(Dur d) {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(d);
}

// Immutable snapshot used by the executor per run
struct CompiledStage {
  // Pointers to manager-owned plan vectors to avoid per-run copies.
  // Lifetime: valid for the duration of a single run(). The manager creates
  // this snapshot after compile(). invalidate_compilation() clears/rebuilds
  // the underlying containers before any subsequent run().
  const std::vector<std::string>* targets = nullptr;
  const std::vector<std::function<std::unique_ptr<Computation<compute::PipelineContext, Mut::ReadWrite>>()>>* factories = nullptr;
  // Pre-resolved labels from targets. Each label stores a string_view that
  // points into the strings in *targets. Because *targets is stable for the
  // duration of a run(), these views remain valid while executing the run.
  std::vector<ComputationLabel> labels;
  bool all_thread_safe = true;
};

template <typename MetricsPolicy>
class StageExecutor {
public:
  using Metrics      = MetricsPolicy;
  using Clock        = typename Metrics::Clock;
  using ThreadHandle = typename Metrics::ThreadHandle;
  using EngineT = ComputeEngine<compute::PipelineContext, Mut::ReadWrite>;

  using PlanBuildHook = void(*)(std::size_t stage_index, std::size_t run_token);

  struct StagePlanCacheEntry {
    ExecOrder plan{};
    bool ready = false;
  };

  struct StageCacheState {
    std::size_t run_seen = 0;
    std::vector<StagePlanCacheEntry> entries;

    void reset() {
      run_seen = 0;
      entries.clear();
    }
  };

  struct ThreadLocalState {
    std::unique_ptr<EngineT> engine;
    compute::PipelineContext data;
    std::size_t run_seen = 0;
    ThreadHandle metrics_handle;
    StageCacheState stage_cache;

    ~ThreadLocalState() {
      StageExecutor::cleanup_hook_(*this);
    }
  };

  using TlsCleanupHook = void(*)(ThreadLocalState&);

  StageExecutor(const CompiledStage& stage,
                ChannelMultiplexer& mux,
                std::size_t run_token,
                Metrics& metrics,
                pipeline::DataFieldSet requirements,
                std::shared_ptr<IRunObserver> observer = {})
    : stage_(stage), mux_(mux), run_token_(run_token), metrics_(metrics), observer_(std::move(observer)), requirements_(requirements) {}

  static void set_tls_cleanup_hook(TlsCleanupHook hook) {
    cleanup_hook_ = hook ? hook : &StageExecutor::default_tls_cleanup;
  }

  static void reset_tls_cleanup_hook() {
    cleanup_hook_ = &StageExecutor::default_tls_cleanup;
  }

  static void set_plan_build_hook(PlanBuildHook hook) {
    plan_build_hook_.store(hook, std::memory_order_relaxed);
  }

  static void reset_plan_build_hook() {
    plan_build_hook_.store(nullptr, std::memory_order_relaxed);
  }

  static void clear_thread_local_state(ThreadLocalState& state) { clear_tls_state(state); }

  template <typename Source>
  void run(Source& src, std::size_t threads) {

    auto logger = Logger::get_logger();
    logger->debug("StageExecutor[run_token={}]: starting run with {} thread(s) and {} stage(s)",
                  run_token_, threads, stage_.labels.size());

    metrics_.configure_stage_breakdown(stage_.labels.size());

    struct InflightGuard {
      Metrics& metrics;
      bool active{false};
      explicit InflightGuard(Metrics& m) : metrics(m), active(true) {
        metrics.on_item_inflight_enter();
      }
      void release() {
        if (!active) return;
        metrics.on_item_inflight_exit();
        active = false;
      }
      ~InflightGuard() { release(); }
    };

    // Worker stage: process each PipelineItem end-to-end on a worker thread
    auto observer = observer_;
    Stage<PipelineItem, void> process_stage(
      [this, observer](PipelineItem item, IEmitter<void>&) {
        auto& tls_state = get_tls_state();
        auto& metrics_handle = ensure_metrics_handle(tls_state);
        InflightGuard inflight_guard{metrics_};

        auto logger = Logger::get_logger();
        const char* session_label = item.session_id.empty() ? "<none>" : item.session_id.c_str();
        logger->debug(
          "StageExecutor[run_token={}]: processing item session='{}' conformer={} path='{}'",
          run_token_, session_label, item.conformer_id, item.item_path);

        metrics_.inc_items_total(metrics_handle);
        if (observer) observer->on_item_begin(run_token_, item);
        const auto stage_begin = Clock::now();

        // Rebuild the worker-local engine once per run() based on the run token.
        // Clear the engine first to release any held resources before rebuilding.
        if (!tls_state.engine || tls_state.run_seen != run_token_) {
          tls_state.engine.reset();  // Explicitly destroy old engine before creating new one
          tls_state.engine = std::make_unique<EngineT>(tls_state.data);
          for (const auto& make : *stage_.factories) {
            tls_state.engine->add(make());
          }
          tls_state.run_seen = run_token_;
        }

        TaskContext ctx;
        tls_state.engine->reset();

        const auto after_setup = Clock::now();
        metrics_.add_setup(metrics_handle, to_ns(after_setup - stage_begin));

        // Per-item temporaries, lifetimes are scoped to this item, and retained by ctx
        std::unique_ptr<StreamSession::Permit> permit;
        std::shared_ptr<RDGeom::POINT3D_VECT>  shared_coords;
        std::shared_ptr<RDKit::Conformer>      conformer;

        const auto prepare_start = after_setup;
        // Prepare per-item state and publish well-known objects into TaskContext
        const bool prepared = prepare_item_state(item, ctx, tls_state.data, permit, shared_coords, conformer, metrics_handle);
        const auto after_prepare = Clock::now();
        metrics_.add_prepare(metrics_handle, to_ns(after_prepare - prepare_start));
        if (!prepared) {
          logger->warn(
            "StageExecutor[run_token={}]: skipping item session='{}' conformer={} due to coordinate/setup failure",
            run_token_, session_label, item.conformer_id);
          metrics_.inc_items_skipped(metrics_handle);
          if (observer) observer->on_item_skipped(run_token_, item, "prepare_failed");
          if (observer) observer->on_item_end(run_token_, item);
          return;
        }

        const auto compute_start = after_prepare;
        for (std::size_t index = 0; index < stage_.labels.size(); ++index) {
          const auto& lbl = stage_.labels[index];
          const auto stage_run_begin = Clock::now();
          auto& plan_entry = ensure_plan_entry(tls_state, index);
          if (!plan_entry.ready) {
            plan_entry.plan = tls_state.engine->plan_from(lbl);
            plan_entry.ready = true;
            notify_plan_built(index);
          }
          tls_state.engine->execute_plan(plan_entry.plan);
          const auto stage_run_end = Clock::now();
          metrics_.add_stage_compute(metrics_handle, index, to_ns(stage_run_end - stage_run_begin));
          auto res = tls_state.engine->get_computation_result(lbl);
          const bool ok = !res.has_error();
          if (observer) observer->on_stage_complete(run_token_, item, lbl.to_string_view(), ok);
          if (!ok) {
            auto& err = res.error();
            logger->warn(
              "StageExecutor[run_token={}]: stage '{}' failed for session='{}' conformer={} path='{}': {}. "
              "Skipping remaining stages for this item.",
              run_token_, lbl.to_string_view(), session_label, item.conformer_id, item.item_path, err.get_message());
            metrics_.inc_items_skipped(metrics_handle);
            if (observer) observer->on_item_skipped(run_token_, item, "stage_failed");
            const auto stage_end = Clock::now();
            metrics_.add_stage_setup(metrics_handle, index, to_ns(stage_end - stage_run_end));
            break;
          }
          if (res.has_value() && res.get_type() == typeid(EmissionList)) {
            auto emits = res.template move_value<EmissionList>();
            for (auto &e : emits) mux_.emit(std::move(e));
          }
          const auto stage_end = Clock::now();
          metrics_.add_stage_setup(metrics_handle, index, to_ns(stage_end - stage_run_end));
        }
        const auto after_compute = Clock::now();
        metrics_.add_compute(metrics_handle, to_ns(after_compute - compute_start));
        if (observer) observer->on_item_end(run_token_, item);
      },
      /*thread_safe=*/stage_.all_thread_safe
    );

    struct Done : IEmitter<void> { void emit() override {} } done;

    struct TimedSource {
      Source& inner;
      StageExecutor& executor;
      Metrics& metrics;

      auto next() {
        const auto start = Clock::now();
        auto item = inner.next();
        const auto end = Clock::now();
        auto& tls_state = StageExecutor::get_tls_state();
        auto& handle = executor.ensure_metrics_handle(tls_state);
        metrics.add_ingest(handle, to_ns(end - start));
        return item;
      }
    } timed_src{src, *this, metrics_};

    PipelineEngine{threads}.run(timed_src, process_stage, done);

    //
    // TLS cleanup: always clear the caller-thread cache so any held resources (including
    // Python objects) release while the embedding runtime is still alive. The registered
    // hook may acquire the GIL when Lahuta is driven from Python. - Besian, October 2025
    //
    cleanup_hook_(get_tls_state());
    logger->debug("StageExecutor[run_token={}]: finished run", run_token_);
  }

private:
  friend struct ThreadLocalState;

  static ThreadLocalState& get_tls_state() {
    static thread_local ThreadLocalState state;
    return state;
  }

  static void clear_tls_state(ThreadLocalState& state) {
    state.metrics_handle.reset();
    state.engine.reset();
    state.data = compute::PipelineContext{};
    state.run_seen = 0;
    state.stage_cache.reset();
  }

  static void default_tls_cleanup(ThreadLocalState& state) {
    clear_tls_state(state);
  }

  void publish_required_handles(dynamic::TaskContext& ctx, compute::PipelineContext& data) const {
    if (!data.session) return;
    if (requirements_.empty()) return;

    auto slices = data.session->model_payload(requirements_);
    if (slices.empty()) return;

    //
    // LMDB-backed views in ModelPayloadSlices pin a read txn on the creating thread.
    // Tasks must copy view data they intend to keep after the task returns or across threads, otherwise
    // dangling views can crash when the underlying txn is released on another thread.
    //
    auto payload = std::make_shared<pipeline::ModelPayloadSlices>(std::move(slices));
    ctx.set_object<pipeline::ModelPayloadSlices>(pipeline::CTX_MODEL_PAYLOAD_KEY, payload);
  }

  ThreadHandle& ensure_metrics_handle(ThreadLocalState& state) {
    metrics_.ensure(state.metrics_handle);
    return state.metrics_handle;
  }

  void notify_plan_built(std::size_t stage_index) const {
    auto hook = plan_build_hook_.load(std::memory_order_relaxed);
    if (hook) hook(stage_index, run_token_);
  }

  StagePlanCacheEntry& ensure_plan_entry(ThreadLocalState& state, std::size_t index) {
    auto& cache = state.stage_cache;
    if (cache.run_seen != run_token_) {
      cache.reset();
      cache.run_seen = run_token_;
    }
    if (cache.entries.size() < stage_.labels.size()) {
      cache.entries.resize(stage_.labels.size());
    }
    return cache.entries[index];
  }

  static inline TlsCleanupHook cleanup_hook_ = &StageExecutor::default_tls_cleanup;
  static inline std::atomic<PlanBuildHook> plan_build_hook_{nullptr};

  //
  // Prepares per-item state and injects it into the task context.
  // Responsibilities:
  // - Bind per-item metadata into the per-thread PipelineContext (path, ids, frame handles)
  // - Write FrameMetadata into TaskContext under CTX_FRAME_KEY
  // - Acquire a StreamSession::Permit to bound in-flight frames (if session present)
  // - Load frame coordinates and construct a lightweight RDKit::Conformer bound to them
  // - Publish coordinates and conformer into TaskContext under CTX_COORDINATES_KEY and CTX_CONFORMER_KEY
  // Returns false if coordinates cannot be prepared. Calling code should skip processing this item.
  //
  bool prepare_item_state(
      const PipelineItem& item,
      TaskContext& t_ctx,
      compute::PipelineContext& p_ctx,
      std::unique_ptr<StreamSession::Permit>& permit,
      std::shared_ptr<RDGeom::POINT3D_VECT>& shared_coords,
      std::shared_ptr<RDKit::Conformer>& conformer,
      ThreadHandle& metrics_handle) {

    p_ctx.ctx          = &t_ctx;
    p_ctx.item_path    = item.item_path;
    p_ctx.conformer_id = item.conformer_id;
    p_ctx.session      = item.session;
    p_ctx.frame        = item.frame;

    compute::set_frame_metadata(t_ctx, item);
    publish_required_handles(t_ctx, p_ctx);

    if (item.session) {
      const auto permit_begin = Clock::now();
      auto acquired = item.session->acquire_permit();
      const auto permit_end = Clock::now();
      permit = std::make_unique<StreamSession::Permit>(std::move(acquired));
      metrics_.add_permit_wait(metrics_handle, to_ns(permit_end - permit_begin));
    } else {
      permit.reset();
    }

    shared_coords.reset();
    conformer.reset();
    if (item.frame) {
      auto view = item.frame->load_coordinates();

      auto positions = view.shared_positions;
      if (!positions) {
        shared_coords = std::make_shared<RDGeom::POINT3D_VECT>(std::move(view.positions));
        if (!shared_coords) return false;
        positions = shared_coords;
      }

      conformer = std::make_shared<RDKit::Conformer>();
      conformer->set3D(true);
      conformer->setId(static_cast<unsigned int>(item.conformer_id));
      conformer->bindExternalPositions(positions);
    }

    if (conformer) {
      t_ctx.set_object<RDKit::Conformer>(pipeline::CTX_CONFORMER_KEY, conformer);
    }
    return true;
  }

  const CompiledStage &stage_;
  ChannelMultiplexer  &mux_;
  std::size_t run_token_;
  Metrics& metrics_;
  std::shared_ptr<IRunObserver> observer_;
  pipeline::DataFieldSet requirements_ = pipeline::DataFieldSet::none();
};

} // namespace lahuta::pipeline::dynamic

#endif // LAHUTA_PIPELINE_DYNAMIC_EXECUTOR_HPP
