#include <atomic>
#include <condition_variable>
#include <mutex>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "pipeline/dynamic/manager.hpp"
#include "pipeline/dynamic/sources.hpp"
#include "pipeline/dynamic/types.hpp"

using namespace lahuta::pipeline::dynamic;

namespace {

// Barrier and concurrency probe shared between worker threads
struct Probe {
  std::mutex m;
  std::condition_variable cv;
  int target = 1;         // number of arrivals to release barrier
  int arrived = 0;        // number of threads that reached barrier
  int current = 0;        // current in-flight in run()
  int max_concurrent = 0; // peak concurrency observed
};

class BarrierProbeTask : public ITask {
public:
  explicit BarrierProbeTask(std::shared_ptr<Probe> p) : p_(std::move(p)) {}

  TaskResult run(const std::string & /*item_path*/, TaskContext & /*ctx*/) override {
    std::unique_lock<std::mutex> lk(p_->m);
    // update concurrency
    p_->current += 1;
    if (p_->current > p_->max_concurrent) p_->max_concurrent = p_->current;

    // wait until target arrivals or bypass if target==1
    p_->arrived += 1;
    if (p_->arrived == p_->target) {
      p_->cv.notify_all();
    } else if (p_->target > 1) {
      p_->cv.wait(lk, [&] { return p_->arrived >= p_->target; });
    }

    p_->current -= 1;
    return TaskResult{true, {}};
  }

private:
  std::shared_ptr<Probe> p_;
};

std::vector<std::string> make_items(int n) {
  std::vector<std::string> v;
  v.reserve(n);
  for (int i = 0; i < n; ++i)
    v.emplace_back("item_" + std::to_string(i));
  return v;
}

} // namespace

TEST(DynamicPipelineParallelism, RunsInParallelAcrossItems) {
  const int threads = 4;
  const int items = threads;

  // Source with N items
  auto src = sources_factory::from_vector(make_items(items));
  StageManager mgr(std::move(src));
  mgr.set_auto_builtins(false); // keep graph minimal

  // Task: thread-safe, uses barrier to ensure all threads arrive before proceeding
  auto probe = std::make_shared<Probe>();
  probe->target = threads;
  auto task = std::make_shared<BarrierProbeTask>(probe);
  mgr.add_task("probe", /*deps=*/{}, task, /*thread_safe=*/true);

  mgr.compile();
  auto report = mgr.run(threads);

  EXPECT_EQ(probe->arrived, items);
  EXPECT_EQ(probe->max_concurrent, threads);
  EXPECT_EQ(report.items_total, static_cast<std::size_t>(items));
  EXPECT_EQ(report.items_processed, static_cast<std::size_t>(items));
  EXPECT_EQ(report.items_skipped, std::size_t{0});
  EXPECT_EQ(report.threads_used, static_cast<std::size_t>(threads));
}

TEST(DynamicPipelineParallelism, CollapsesWhenTaskMarkedUnsafe) {
  const int threads = 4;
  const int items = threads * 2;

  auto src2 = sources_factory::from_vector(make_items(items));
  StageManager mgr(std::move(src2));
  mgr.set_auto_builtins(false);

  auto probe = std::make_shared<Probe>();
  probe->target = 1; // do not wait for others, expect serialized execution
  auto task = std::make_shared<BarrierProbeTask>(probe);
  mgr.add_task("probe", /*deps=*/{}, task, /*thread_safe=*/false);

  mgr.compile();
  auto report = mgr.run(threads); // ask for 4, but expect single-threaded due to unsafe task

  EXPECT_EQ(probe->arrived, items);
  EXPECT_EQ(probe->max_concurrent, 1);
  EXPECT_EQ(report.threads_used, std::size_t{1});
}

namespace {
using namespace lahuta::topology::compute;
using namespace lahuta::pipeline::compute;

struct TestParams : ParameterBase<TestParams> {
  static constexpr ParameterInterface::TypeId TYPE_ID = 250;
};

// Simple counting computation. Increments a shared counter when executed.
class CountingComputation : public Computation<PipelineContext, Mut::ReadWrite> {
public:
  CountingComputation(
      std::string label, std::shared_ptr<std::atomic<int>> counter, std::vector<ComputationLabel> deps = {})
      : label_storage_(std::move(label)), label_(label_storage_), counter_(std::move(counter)),
        deps_(std::move(deps)) {}

  ComputationResult
  execute(DataContext<PipelineContext, Mut::ReadWrite> &, const ParameterInterface &p) override {
    if (p.type_id() != TestParams::TYPE_ID) return ComputationResult(ComputationError("bad params"));
    ++(*counter_);
    return ComputationResult(true);
  }
  std::unique_ptr<ParameterInterface> get_parameters() const override {
    return std::make_unique<TestParams>();
  }
  const ComputationLabel &get_label() const override { return label_; }
  std::vector<ComputationLabel> get_dependencies() const override { return deps_; }

private:
  std::string label_storage_;
  ComputationLabel label_;
  std::shared_ptr<std::atomic<int>> counter_;
  std::vector<ComputationLabel> deps_;
};
} // namespace

TEST(DynamicPipelineParallelism, MemoizationAcrossTargetsRunsOncePerItem) {
  const int threads = 3;
  const int items = 7;

  auto src3 = sources_factory::from_vector(make_items(items));
  StageManager mgr(std::move(src3));
  mgr.set_auto_builtins(false);

  auto cntA = std::make_shared<std::atomic<int>>(0);
  auto cntB = std::make_shared<std::atomic<int>>(0);

  mgr.add_computation(
      "A",
      /*deps=*/{},
      [cntA]() { return std::make_unique<CountingComputation>("A", cntA); },
      /*thread_safe=*/true);

  mgr.add_computation(
      "B",
      /*deps=*/{"A"},
      [cntB]() {
        std::vector<ComputationLabel> d{ComputationLabel{"A"}};
        return std::make_unique<CountingComputation>("B", cntB, d);
      },
      /*thread_safe=*/true);

  mgr.compile();
  mgr.run(threads);

  // A is a dependency of B and also a target. run_from should execute A once per item (memoized) not twice.
  EXPECT_EQ(cntA->load(), items);
  EXPECT_EQ(cntB->load(), items);
}
