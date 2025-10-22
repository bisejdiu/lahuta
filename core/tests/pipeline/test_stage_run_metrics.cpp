#include <chrono>
#include <cstdint>
#include <thread>
#include <vector>

#include <gtest/gtest.h>

#include "pipeline/dynamic/run_metrics.hpp"

using namespace lahuta::pipeline::dynamic;

TEST(StageRunMetricsTest, AggregatesSingleThread) {
  StageRunMetrics metrics;
  StageRunMetrics::ThreadHandle handle;
  metrics.ensure(handle);

  metrics.add_ingest(handle, std::chrono::nanoseconds(5));
  metrics.add_prepare(handle, std::chrono::nanoseconds(3));
  metrics.add_setup(handle, std::chrono::nanoseconds(7));
  metrics.add_compute(handle, std::chrono::nanoseconds(11));
  metrics.inc_items_total(handle);
  metrics.inc_items_skipped(handle);

  metrics.add_flush(std::chrono::nanoseconds(13));
  metrics.on_item_inflight_enter();
  metrics.on_item_inflight_exit();

  const auto snapshot = metrics.snapshot();
  EXPECT_EQ(snapshot.ingest_ns, 5);
  EXPECT_EQ(snapshot.prepare_ns, 3);
  EXPECT_EQ(snapshot.setup_ns, 7);
  EXPECT_EQ(snapshot.compute_ns, 11);
  EXPECT_EQ(snapshot.flush_ns, 13);
  EXPECT_EQ(snapshot.items_total, std::size_t{1});
  EXPECT_EQ(snapshot.items_skipped, std::size_t{1});
  EXPECT_EQ(snapshot.inflight_peak, std::size_t{1});
  EXPECT_EQ(snapshot.inflight_samples, std::uint64_t{2});
  EXPECT_EQ(snapshot.inflight_sum, std::uint64_t{1});
  EXPECT_EQ(snapshot.permit_wait_ns_total, std::uint64_t{0});
}

TEST(StageRunMetricsTest, AggregatesAcrossThreads) {
  StageRunMetrics metrics;

  constexpr int threads = 4;

  std::vector<std::thread> workers;
  workers.reserve(threads);

  for (int i = 0; i < threads; ++i) {
    workers.emplace_back([&metrics, i]() {
      StageRunMetrics::ThreadHandle handle;
      metrics.ensure(handle);
      metrics.add_ingest(handle, std::chrono::nanoseconds(10 + i));
      metrics.add_prepare(handle, std::chrono::nanoseconds(20 + i));
      metrics.add_setup(handle, std::chrono::nanoseconds(30 + i));
      metrics.add_compute(handle, std::chrono::nanoseconds(40 + i));
      metrics.inc_items_total(handle);
      if ((i % 2) == 0) {
        metrics.inc_items_skipped(handle);
      }
    });
  }

  for (auto& t : workers) {
    t.join();
  }

  metrics.add_flush(std::chrono::nanoseconds(100));
  metrics.on_item_inflight_enter();
  metrics.on_item_inflight_enter();
  metrics.on_item_inflight_exit();
  metrics.on_item_inflight_exit();

  std::int64_t expected_ingest = 0;
  std::int64_t expected_prepare = 0;
  std::int64_t expected_setup = 0;
  std::int64_t expected_compute = 0;
  std::size_t expected_skipped = 0;
  for (int i = 0; i < threads; ++i) {
    expected_ingest  += 10 + i;
    expected_prepare += 20 + i;
    expected_setup   += 30 + i;
    expected_compute += 40 + i;
    if ((i % 2) == 0) {
      ++expected_skipped;
    }
  }

  const auto snapshot = metrics.snapshot();
  EXPECT_EQ(snapshot.ingest_ns, expected_ingest);
  EXPECT_EQ(snapshot.prepare_ns, expected_prepare);
  EXPECT_EQ(snapshot.setup_ns, expected_setup);
  EXPECT_EQ(snapshot.compute_ns, expected_compute);
  EXPECT_EQ(snapshot.flush_ns, 100);
  EXPECT_EQ(snapshot.items_total, static_cast<std::size_t>(threads));
  EXPECT_EQ(snapshot.items_skipped, expected_skipped);
  EXPECT_EQ(snapshot.inflight_peak, std::size_t{2});
  EXPECT_EQ(snapshot.permit_wait_ns_total, std::uint64_t{0});
}

TEST(StageRunMetricsTest, RecordsStageBreakdownWhenEnabled) {
  StageRunMetrics metrics(/*enable_stage_breakdown=*/true);
  metrics.configure_stage_breakdown(2);
  StageRunMetrics::ThreadHandle handle;
  metrics.ensure(handle);
  metrics.add_stage_setup(handle, 0, std::chrono::nanoseconds(10));
  metrics.add_stage_compute(handle, 0, std::chrono::nanoseconds(20));
  metrics.add_stage_setup(handle, 1, std::chrono::nanoseconds(30));
  metrics.add_stage_compute(handle, 1, std::chrono::nanoseconds(40));
  metrics.add_permit_wait(handle, std::chrono::nanoseconds(50));
  metrics.add_permit_wait(handle, std::chrono::nanoseconds(150));

  const auto snapshot = metrics.snapshot();
  ASSERT_EQ(snapshot.stage_setup_ns.size(), std::size_t{2});
  ASSERT_EQ(snapshot.stage_compute_ns.size(), std::size_t{2});
  EXPECT_EQ(snapshot.stage_setup_ns[0], 10);
  EXPECT_EQ(snapshot.stage_compute_ns[0], 20);
  EXPECT_EQ(snapshot.stage_setup_ns[1], 30);
  EXPECT_EQ(snapshot.stage_compute_ns[1], 40);
  EXPECT_EQ(snapshot.permit_wait_ns_total, std::uint64_t{200});
  EXPECT_EQ(snapshot.permit_wait_samples, std::uint64_t{2});
  EXPECT_EQ(snapshot.permit_wait_ns_min, std::int64_t{50});
  EXPECT_EQ(snapshot.permit_wait_ns_max, std::int64_t{150});
}
