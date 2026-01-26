#include <atomic>
#include <optional>
#include <vector>

#include <gtest/gtest.h>

#include "compute/compute_base.hpp"
#include "compute/parameters.hpp"
#include "pipeline/data/pipeline_item.hpp"
#include "pipeline/io/channel_multiplexer.hpp"
#include "pipeline/metrics/run_metrics.hpp"
#include "pipeline/runtime/executor.hpp"
#include "pipeline/task/compute/context.hpp"

using namespace lahuta;
namespace P = lahuta::pipeline;

namespace {

struct CountingParams : public P::ParameterBase<CountingParams> {
  static constexpr P::ParameterInterface::TypeId TYPE_ID = 246;
};

class CountingComputation : public P::Computation<P::PipelineContext, P::Mut::ReadWrite> {
public:
  explicit CountingComputation(std::string label) : label_buffer_(std::move(label)), label_(label_buffer_) {}

  std::unique_ptr<P::ParameterInterface> get_parameters() const override {
    return std::make_unique<CountingParams>();
  }

  const P::ComputationLabel &get_label() const override { return label_; }
  std::vector<P::ComputationLabel> get_dependencies() const override { return {}; }

  P::ComputationResult execute(P::DataContext<P::PipelineContext, P::Mut::ReadWrite> &ctx,
                               const P::ParameterInterface &) override {
    (void)ctx;
    exec_calls_.fetch_add(1, std::memory_order_relaxed);
    return P::ComputationResult(true);
  }

  static void reset() { exec_calls_.store(0, std::memory_order_relaxed); }
  static int run_count() { return exec_calls_.load(std::memory_order_relaxed); }

private:
  std::string label_buffer_;
  P::ComputationLabel label_;
  inline static std::atomic<int> exec_calls_{0};
};

struct VectorSource {
  std::vector<P::PipelineItem> items;
  std::size_t index = 0;

  std::optional<P::PipelineItem> next() {
    if (index >= items.size()) return std::nullopt;
    return items[index++];
  }
};

std::atomic<int> g_plan_builds{0};

void record_plan_build(std::size_t, std::size_t) { g_plan_builds.fetch_add(1, std::memory_order_relaxed); }

VectorSource make_source(std::size_t count) {
  VectorSource src;
  src.items.resize(count);
  src.index = 0;
  return src;
}

} // namespace

TEST(StageExecutorPlanCacheTest, BuildsEachPlanExactlyOncePerRunToken) {
  using Executor = P::StageExecutor<P::NullStageRunMetrics>;

  std::vector<std::string> targets = {"count_stage"};
  std::vector<std::function<std::unique_ptr<P::Computation<P::PipelineContext, P::Mut::ReadWrite>>()>>
      factories;
  factories.emplace_back([] { return std::make_unique<CountingComputation>("count_stage"); });

  P::CompiledStage stage{};
  stage.targets         = &targets;
  stage.factories       = &factories;
  stage.labels          = {P::ComputationLabel{targets[0].c_str()}};
  stage.all_thread_safe = true;

  P::ChannelMultiplexer mux;
  P::NullStageRunMetrics metrics;
  const auto requirements = P::DataFieldSet::none();

  g_plan_builds.store(0, std::memory_order_relaxed);
  CountingComputation::reset();
  Executor::set_plan_build_hook(&record_plan_build);

  {
    VectorSource src = make_source(5);
    Executor executor(stage, mux, /*run_token=*/1, metrics, requirements);
    executor.run(src, /*threads=*/1);
  }

  EXPECT_EQ(g_plan_builds.load(std::memory_order_relaxed), 1);
  EXPECT_EQ(CountingComputation::run_count(), 5);

  {
    VectorSource src = make_source(3);
    Executor executor(stage, mux, /*run_token=*/2, metrics, requirements);
    executor.run(src, /*threads=*/1);
  }

  EXPECT_EQ(g_plan_builds.load(std::memory_order_relaxed), 2);
  EXPECT_EQ(CountingComputation::run_count(), 8);

  Executor::reset_plan_build_hook();
}
