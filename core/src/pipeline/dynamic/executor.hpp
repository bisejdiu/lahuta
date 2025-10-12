#ifndef LAHUTA_PIPELINE_DYNAMIC_EXECUTOR_HPP
#define LAHUTA_PIPELINE_DYNAMIC_EXECUTOR_HPP

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
#include "pipeline/dynamic/channel_multiplexer.hpp"
#include "pipeline/dynamic/keys.hpp"
#include "pipeline/dynamic/types.hpp"
#include "pipeline/engine.hpp"
#include "pipeline/pipeline_item.hpp"
#include "pipeline/stream_session.hpp"

// clang-format off
namespace lahuta::pipeline::dynamic {
using namespace lahuta::topology::compute;

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

// Per-run executor. Reuse a thread-local ComputeEngine within a single run() and resets it per item.
class StageExecutor {
public:
  using EngineT = ComputeEngine<compute::PipelineContext, Mut::ReadWrite>;

  struct ThreadLocalState {
    std::unique_ptr<EngineT> engine;
    compute::PipelineContext data;
    std::size_t run_seen = 0;

    ~ThreadLocalState() {
      StageExecutor::cleanup_hook_(*this);
    }
  };

  using TlsCleanupHook = void(*)(ThreadLocalState&);

  StageExecutor(const CompiledStage& stage, ChannelMultiplexer& mux, std::size_t run_token)
    : stage_(stage), mux_(mux), run_token_(run_token) {}

  static void set_tls_cleanup_hook(TlsCleanupHook hook) {
    cleanup_hook_ = hook ? hook : &StageExecutor::default_tls_cleanup;
  }

  static void reset_tls_cleanup_hook() {
    cleanup_hook_ = &StageExecutor::default_tls_cleanup;
  }

  static void clear_thread_local_state(ThreadLocalState& state) {
    clear_tls_state(state);
  }

  template <typename Source>
  void run(Source& src, std::size_t threads) {

    auto logger = Logger::get_logger();
    logger->debug("StageExecutor[run_token={}]: starting run with {} thread(s) and {} stage(s)",
                  run_token_, threads, stage_.labels.size());

    // Worker stage: process each PipelineItem end-to-end on a worker thread
    Stage<PipelineItem, void> process_stage(
      [this](PipelineItem item, IEmitter<void>&) {
        auto& tls_state = get_tls_state();

        auto logger = Logger::get_logger();
        const char* session_label = item.session_id.empty() ? "<none>" : item.session_id.c_str();
        logger->debug(
          "StageExecutor[run_token={}]: processing item session='{}' conformer={} path='{}'",
          run_token_, session_label, item.conformer_id, item.item_path);

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

        // Per-item temporaries, lifetimes are scoped to this item, and retained by ctx
        std::unique_ptr<StreamSession::Permit> permit;
        std::shared_ptr<RDGeom::POINT3D_VECT>  shared_coords;
        std::shared_ptr<RDKit::Conformer>      conformer;

        // Prepare per-item state and publish well-known objects into TaskContext
        if (!prepare_item_state(item, ctx, tls_state.data, permit, shared_coords, conformer)) {
          logger->warn(
            "StageExecutor[run_token={}]: skipping item session='{}' conformer={} due to coordinate/setup failure",
            run_token_, session_label, item.conformer_id);
          return;
        }

        for (const auto& lbl : stage_.labels) {
          (void)tls_state.engine->template run_from<void>(lbl);
          auto res = tls_state.engine->get_computation_result(lbl);
          if (res.has_error()) break;
          if (res.has_value() && res.get_type() == typeid(EmissionList)) {
            auto emits = res.template move_value<EmissionList>();
            for (auto &e : emits) mux_.emit(std::move(e));
          }
        }
      },
      /*thread_safe=*/stage_.all_thread_safe
    );

    struct Done : IEmitter<void> { void emit() override {} } done;
    PipelineEngine{threads}.run(src, process_stage, done);

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
    state.engine.reset();
    state.data = compute::PipelineContext{};
    state.run_seen = 0;
  }

  static void default_tls_cleanup(ThreadLocalState& state) {
    clear_tls_state(state);
  }

  static inline TlsCleanupHook cleanup_hook_ = &StageExecutor::default_tls_cleanup;

  // Prepares per-item state and injects it into the task context.
  // Responsibilities:
  // - Bind per-item metadata into the per-thread PipelineContext (path, ids, frame handles)
  // - Write FrameMetadata into TaskContext under CTX_FRAME_KEY
  // - Acquire a StreamSession::Permit to bound in-flight frames (if session present)
  // - Load frame coordinates and construct a lightweight RDKit::Conformer bound to them
  // - Publish coordinates and conformer into TaskContext under CTX_COORDINATES_KEY and CTX_CONFORMER_KEY
  // Returns false if coordinates cannot be prepared. Calling code should skip processing this item.
  bool prepare_item_state(
      const PipelineItem& item,
      TaskContext& t_ctx,
      compute::PipelineContext& p_ctx,
      std::unique_ptr<StreamSession::Permit>& permit,
      std::shared_ptr<RDGeom::POINT3D_VECT>& shared_coords,
      std::shared_ptr<RDKit::Conformer>& conformer) {

    p_ctx.ctx          = &t_ctx;
    p_ctx.item_path    = item.item_path;
    p_ctx.conformer_id = item.conformer_id;
    p_ctx.session      = item.session;
    p_ctx.frame        = item.frame;

    compute::set_frame_metadata(t_ctx, item);

    if (item.session) {
      permit = std::make_unique<StreamSession::Permit>(item.session->acquire_permit());
    } else {
      permit.reset();
    }

    // Detached conformers
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
};

} // namespace lahuta::pipeline::dynamic

#endif // LAHUTA_PIPELINE_DYNAMIC_EXECUTOR_HPP
