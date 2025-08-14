#ifndef LAHUTA_PIPELINE_DYNAMIC_EXECUTOR_HPP
#define LAHUTA_PIPELINE_DYNAMIC_EXECUTOR_HPP

#include <cstddef>
#include <functional>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "compute/engine.hpp"
#include "pipeline/compute/context.hpp"
#include "pipeline/core/emitter.hpp"
#include "pipeline/core/stage.hpp"
#include "pipeline/dynamic/channel_multiplexer.hpp"
#include "pipeline/dynamic/types.hpp"
#include "pipeline/engine.hpp"

// clang-format off
namespace lahuta::pipeline::dynamic {
using namespace lahuta::topology::compute;

// The goal of CompiledStage is to avoid per-run copies
// Immutable snapshot used by the executor per run
struct CompiledStage {
  // Pointers to manager-owned plan vectors to avoid per-run copies.
  // Lifetime: valid for the duration of a single run(). The manager creates
  // this snapshot after compile(). invalidate_compilation() clears/rebuilds
  // the underlying containers before any subsequent run() ----> so no stale access (hopefully?)
  const std::vector<std::string>* targets = nullptr;
  const std::vector<std::function<std::unique_ptr<Computation<compute::PipelineContext, Mut::ReadWrite>>()>>* factories = nullptr;
  // Pre-resolved labels from targets. Each label stores a string_view that
  // points into the strings in *targets. Because *targets is stable for the
  // duration of a run(), these views remain valid while executing the run.
  std::vector<ComputationLabel> labels;
  bool all_thread_safe = true;
};

// Per-run executor. Reuses a thread-local ComputeEngine within a single run() and resets it per item.
class StageExecutor {
public:
  StageExecutor(const CompiledStage& stage, ChannelMultiplexer& mux, std::size_t run_token)
    : stage_(stage), mux_(mux), run_token_(run_token) {}

  template <typename SourceVariant>
  void run(SourceVariant& src, std::size_t threads) {
    using InT = std::string;

    pipeline::Stage<InT, void> process_stage(
      [this](InT path, pipeline::IEmitter<void>&) {
        using EngineT = ComputeEngine<compute::PipelineContext, Mut::ReadWrite>;
        static thread_local std::unique_ptr<EngineT> tls_engine;
        static thread_local compute::PipelineContext tls_data;
        static thread_local std::size_t tls_run_seen = 0;

        // Rebuild the worker-local engine once per run() based on the run token.
        if (!tls_engine || tls_run_seen != run_token_) {
          tls_engine = std::make_unique<EngineT>(tls_data);
          tls_engine->set_auto_heal(true);
          for (const auto& make : *stage_.factories) tls_engine->add(make());
          tls_run_seen = run_token_;
        }

        TaskContext ctx;
        tls_data.item_path = std::move(path);
        tls_data.ctx = &ctx;
        tls_engine->reset();

        for (const auto& lbl : stage_.labels) {
          (void)tls_engine->run_from<void>(lbl);
          auto res = tls_engine->get_computation_result(lbl);
          if (res.has_error()) break;
          if (res.has_value() && res.get_type() == typeid(dynamic::EmissionList)) {
            auto emits = res.move_value<dynamic::EmissionList>();
            for (auto &e : emits) mux_.emit(std::move(e));
          }
        }
      },
      /*thread_safe=*/stage_.all_thread_safe
    );

    struct Done : pipeline::IEmitter<void> { void emit() override {} } done;

    std::visit([&](auto& s) {
      pipeline::PipelineEngine{threads}.run(s, process_stage, done);
    }, src);
  }

private:
  const CompiledStage& stage_;
  ChannelMultiplexer& mux_;
  std::size_t run_token_;
};

} // namespace lahuta::pipeline::dynamic

#endif // LAHUTA_PIPELINE_DYNAMIC_EXECUTOR_HPP
