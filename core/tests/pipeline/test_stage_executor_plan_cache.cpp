#include <atomic>
#include <optional>
#include <vector>

#include <gtest/gtest.h>

#include "compute/compute_base.hpp"
#include "compute/parameters.hpp"
#include "pipeline/compute/context.hpp"
#include "pipeline/dynamic/channel_multiplexer.hpp"
#include "pipeline/dynamic/executor.hpp"
#include "pipeline/dynamic/run_metrics.hpp"
#include "pipeline/pipeline_item.hpp"

using namespace lahuta;
using namespace lahuta::pipeline::compute;
using namespace lahuta::pipeline::dynamic;

// clang-format off
namespace {

struct CountingParams : public ParameterBase<CountingParams> {
  static constexpr ParameterInterface::TypeId TYPE_ID = 246;
};

class CountingComputation : public Computation<PipelineContext, Mut::ReadWrite> {
public:
  explicit CountingComputation(std::string label) : label_buffer_(std::move(label)), label_(label_buffer_) {}

  std::unique_ptr<ParameterInterface> get_parameters() const override {
    return std::make_unique<CountingParams>();
  }

  const ComputationLabel &get_label() const override { return label_; }
  std::vector<ComputationLabel> get_dependencies() const override { return {}; }

  ComputationResult execute(DataContext<PipelineContext, Mut::ReadWrite> &ctx, const ParameterInterface &) override {
    (void)ctx;
    exec_calls_.fetch_add(1, std::memory_order_relaxed);
    return ComputationResult(true);
  }

  static void reset() { exec_calls_.store(0, std::memory_order_relaxed); }
  static int run_count() { return exec_calls_.load(std::memory_order_relaxed); }

private:
  std::string label_buffer_;
  ComputationLabel label_;
  inline static std::atomic<int> exec_calls_{0};
};

struct VectorSource {
  std::vector<PipelineItem> items;
  std::size_t index = 0;

  std::optional<PipelineItem> next() {
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
  using Executor = StageExecutor<NullStageRunMetrics>;

  std::vector<std::string> targets = {"count_stage"};
  std::vector<std::function<std::unique_ptr<Computation<PipelineContext, Mut::ReadWrite>>()>> factories;
  factories.emplace_back([] { return std::make_unique<CountingComputation>("count_stage"); });

  CompiledStage stage{};
  stage.targets = &targets;
  stage.factories = &factories;
  stage.labels = {ComputationLabel{targets[0].c_str()}};
  stage.all_thread_safe = true;

  ChannelMultiplexer mux;
  NullStageRunMetrics metrics;
  const auto requirements = pipeline::DataFieldSet::none();

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
