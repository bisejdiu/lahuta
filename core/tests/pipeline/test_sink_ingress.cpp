/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s(22, '\0'); char* ptr = s.data();
 *   for (char c : std::string_view{"besian"}) *ptr++ = c;
 *   for (char c : std::string_view{"sejdiu"}) *ptr++ = c;
 *   *ptr++ = '@';
 *   for (char c : std::string_view{"gmail.com"}) *ptr++ = c;
 *   return s;
 * }();
 *
 */

#include <atomic>
#include <chrono>
#include <memory>
#include <mutex>
#include <set>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include <gtest/gtest.h>

#include "pipeline/io/backpressure.hpp"
#include "pipeline/io/sink_iface.hpp"

using namespace lahuta::pipeline;
namespace {

class RecordingSink : public IDynamicSink { // records writes and flags flush/close calls.
public:
  void write(EmissionView v) override {
    std::lock_guard<std::mutex> lk(m_);
    writes.emplace_back(std::string(v.payload));
  }
  void flush() override { flushed.store(true, std::memory_order_release); }
  void close() override { closed.store(true, std::memory_order_release); }

  std::vector<std::string> writes;
  std::atomic<bool> flushed{false};
  std::atomic<bool> closed{false};

private:
  std::mutex m_;
};

BackpressureConfig small_cfg() {
  BackpressureConfig cfg;
  cfg.max_queue_msgs   = 16;
  cfg.max_queue_bytes  = 1 * 1024 * 1024;
  cfg.max_batch_msgs   = 8;
  cfg.max_batch_bytes  = 1 * 1024 * 1024;
  cfg.offer_wait_slice = std::chrono::milliseconds{5};
  cfg.on_full          = OnFull::Block;
  cfg.required         = true;
  return cfg;
}

TEST(SinkIngress, BasicDrainAndJoin) {
  auto sink = std::make_shared<RecordingSink>();
  SinkIngress ingress(sink, small_cfg());
  ingress.start();

  auto payload             = std::make_shared<const std::string>(std::string(128, 'x'));
  constexpr uint32_t ch_id = 1;
  const std::size_t count  = 10;
  for (std::size_t i = 0; i < count; ++i) {
    ASSERT_TRUE(ingress.offer(ch_id, payload, payload->size()));
  }

  // Signal end of stream and join with a reasonable deadline
  ingress.close_queue();
  const auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(5);
  ASSERT_TRUE(ingress.join_until(deadline));

  EXPECT_EQ(sink->writes.size(), count);
  EXPECT_TRUE(sink->flushed.load());
  EXPECT_TRUE(sink->closed.load());
}

TEST(SinkIngress, JoinWithoutOffers) {
  auto sink = std::make_shared<RecordingSink>();
  SinkIngress ingress(sink, small_cfg());
  ingress.start();

  ingress.close_queue();
  const auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(5);
  ASSERT_TRUE(ingress.join_until(deadline));

  EXPECT_TRUE(sink->flushed.load());
  EXPECT_TRUE(sink->closed.load());
  EXPECT_TRUE(sink->writes.empty());
}

TEST(SinkIngress, OfferAfterCloseIsRejected) {
  auto sink = std::make_shared<RecordingSink>();
  SinkIngress ingress(sink, small_cfg());
  ingress.start();

  ingress.close_queue();
  auto payload = std::make_shared<const std::string>("hello");
  // Offering after close should fail
  EXPECT_FALSE(ingress.offer(1, payload, payload->size()));

  const auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(5);
  ASSERT_TRUE(ingress.join_until(deadline));
}

TEST(SinkIngress, RapidJoinAfterClose_NoTimeout) {
  //
  // Regression test for the lost-wakeup race condition fixed in join_until().
  // This test repeatedly performs close() followed immediately by join_until()
  // to trigger the race window where the writer thread might finish and notify
  // before the join thread enters the wait.
  //
  constexpr int ITERATIONS = 100;
  for (int iteration = 0; iteration < ITERATIONS; ++iteration) {
    SCOPED_TRACE(testing::Message() << "Iteration " << iteration);

    auto sink = std::make_shared<RecordingSink>();
    SinkIngress ingress(sink, small_cfg());
    ingress.start();

    // Minimal work - just close and join immediately
    ingress.close_queue();

    //
    // 200ms allows for scheduler delays on busy systems while correct code completes in <5ms.
    //
    const auto deadline = std::chrono::steady_clock::now() + std::chrono::milliseconds(200);
    ASSERT_TRUE(ingress.join_until(deadline)) << "Join timed out (potential race condition reintroduced)";

    EXPECT_TRUE(sink->flushed.load());
    EXPECT_TRUE(sink->closed.load());
  }
}

TEST(SinkIngress, ConcurrentOffersAndClose) {
  // Test that concurrent offers and close don't cause crashes or deadlocks.
  // The key test is that the system handles concurrent access safely.
  auto sink = std::make_shared<RecordingSink>();
  SinkIngress ingress(sink, small_cfg());
  ingress.start();

  std::atomic<int> successful_offers{0};
  std::atomic<int> failed_offers{0};
  constexpr int TotalOffers = 100;

  std::thread producer([&]() {
    auto payload = std::make_shared<const std::string>("data");
    for (int i = 0; i < TotalOffers; ++i) {
      if (ingress.offer(1, payload, payload->size())) {
        ++successful_offers;
      } else {
        ++failed_offers;
      }
      // Small delay to spread offers over time
      std::this_thread::sleep_for(std::chrono::microseconds(10));
    }
  });

  // Let some offers go through before closing
  std::this_thread::sleep_for(std::chrono::milliseconds(50));
  ingress.close_queue();

  producer.join();

  const auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(5);
  ASSERT_TRUE(ingress.join_until(deadline));

  // Verify accounting: all offers were either successful or failed
  EXPECT_EQ(successful_offers.load() + failed_offers.load(), TotalOffers);

  // At least some offers should have succeeded (unless timing is extremely unlucky)
  EXPECT_GT(successful_offers.load(), 0) << "No successful offers";

  // Written count should match successful offers
  EXPECT_EQ(sink->writes.size(), static_cast<std::size_t>(successful_offers.load()));

  // The sink should have been properly closed
  EXPECT_TRUE(sink->flushed.load());
  EXPECT_TRUE(sink->closed.load());
}

TEST(SinkIngress, MultipleJoinsSucceed) {
  // Verify that calling join_until() multiple times is safe (idempotent).
  auto sink = std::make_shared<RecordingSink>();
  SinkIngress ingress(sink, small_cfg());
  ingress.start();

  ingress.close_queue();
  const auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(5);

  // First join
  ASSERT_TRUE(ingress.join_until(deadline));
  EXPECT_TRUE(sink->flushed.load());
  EXPECT_TRUE(sink->closed.load());

  // Second join should also succeed (idempotent behavior)
  const auto deadline2 = std::chrono::steady_clock::now() + std::chrono::seconds(5);
  ASSERT_TRUE(ingress.join_until(deadline2));

  // State should remain consistent
  EXPECT_TRUE(sink->flushed.load());
  EXPECT_TRUE(sink->closed.load());
}

TEST(SinkIngress, JoinWithVeryShortDeadlineTimesOut) {
  // Test that join_until() properly respects the deadline and returns false on timeout.
  auto sink = std::make_shared<RecordingSink>();
  SinkIngress ingress(sink, small_cfg());
  ingress.start();

  // Enqueue some data but don't close the queue - writer will block waiting for more
  auto payload = std::make_shared<const std::string>(std::string(128, 'x'));
  for (int i = 0; i < 5; ++i) {
    ASSERT_TRUE(ingress.offer(1, payload, payload->size()));
  }

  // Try to join with a very short deadline while queue is still open
  const auto deadline = std::chrono::steady_clock::now() + std::chrono::milliseconds(10);
  EXPECT_FALSE(ingress.join_until(deadline)) << "Expected timeout but join succeeded";

  // Now close and join properly
  ingress.close_queue();
  const auto deadline2 = std::chrono::steady_clock::now() + std::chrono::seconds(5);
  ASSERT_TRUE(ingress.join_until(deadline2));
}

class ThreadTrackingSink : public IDynamicSink {
public:
  void write(EmissionView) override {
    {
      std::lock_guard<std::mutex> lk(m_);
      writer_threads_.insert(std::this_thread::get_id());
    }
    write_count_.fetch_add(1, std::memory_order_relaxed);
    std::this_thread::sleep_for(std::chrono::milliseconds(2));
  }

  void flush() override { flushed_.store(true, std::memory_order_release); }
  void close() override { closed_.store(true, std::memory_order_release); }

  std::size_t thread_count() const {
    std::lock_guard<std::mutex> lk(m_);
    return writer_threads_.size();
  }

  std::size_t write_count() const { return write_count_.load(std::memory_order_relaxed); }
  bool flushed() const { return flushed_.load(std::memory_order_acquire); }
  bool closed() const { return closed_.load(std::memory_order_acquire); }

private:
  mutable std::mutex m_;
  std::set<std::thread::id> writer_threads_;
  std::atomic<std::size_t> write_count_{0};
  std::atomic<bool> flushed_{false};
  std::atomic<bool> closed_{false};
};

class FailOnFirstWriteSink : public IDynamicSink {
public:
  void write(EmissionView) override {
    bool already = fail_.exchange(true, std::memory_order_acq_rel);
    if (!already) {
      throw std::runtime_error("Oh, no!");
    }
  }

  void flush() override {}
  void close() override {}

  bool threw() const { return fail_.load(std::memory_order_acquire); }

private:
  std::atomic<bool> fail_{false};
};

TEST(SinkIngress, SingleWriterUsesDefaultThreadCount) {
  auto sink          = std::make_shared<ThreadTrackingSink>();
  auto cfg           = small_cfg();
  cfg.max_batch_msgs = 1;
  SinkIngress ingress(sink, cfg);
  ingress.start();

  auto payload             = std::make_shared<const std::string>("x");
  constexpr int total_msgs = 32;
  for (int i = 0; i < total_msgs; ++i) {
    ASSERT_TRUE(ingress.offer(1, payload, payload->size()));
  }

  ingress.close_queue();
  const auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(5);
  ASSERT_TRUE(ingress.join_until(deadline));

  EXPECT_EQ(sink->write_count(), static_cast<std::size_t>(total_msgs));
  EXPECT_EQ(sink->thread_count(), 1u);
  EXPECT_TRUE(sink->flushed());
  EXPECT_TRUE(sink->closed());
}

TEST(SinkIngress, ConfigurableWriterThreads) {
  auto sink          = std::make_shared<ThreadTrackingSink>();
  auto cfg           = small_cfg();
  cfg.writer_threads = 3;
  cfg.max_batch_msgs = 1;

  SinkIngress ingress(sink, cfg);
  ingress.start();

  constexpr int total_msgs = 90;
  auto payload             = std::make_shared<const std::string>("x");
  for (int i = 0; i < total_msgs; ++i) {
    ASSERT_TRUE(ingress.offer(1, payload, payload->size()));
  }

  ingress.close_queue();
  const auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(5);
  ASSERT_TRUE(ingress.join_until(deadline));

  EXPECT_EQ(sink->write_count(), static_cast<std::size_t>(total_msgs));
  EXPECT_GE(sink->thread_count(), 2u);
  EXPECT_LE(sink->thread_count(), cfg.writer_threads);
  EXPECT_TRUE(sink->flushed());
  EXPECT_TRUE(sink->closed());
}

TEST(SinkIngress, WriterFailureClosesQueue) {
  auto sink          = std::make_shared<FailOnFirstWriteSink>();
  auto cfg           = small_cfg();
  cfg.writer_threads = 4;
  cfg.max_batch_msgs = 1;

  SinkIngress ingress(sink, cfg);
  ingress.start();

  auto payload = std::make_shared<const std::string>("y");
  // Enqueue several messages, the first write will throw
  for (int i = 0; i < 8; ++i) {
    ASSERT_TRUE(ingress.offer(1, payload, payload->size()));
  }

  ingress.close_queue();
  const auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(5);
  ASSERT_TRUE(ingress.join_until(deadline));
  EXPECT_TRUE(sink->threw());

  auto another = std::make_shared<const std::string>("z");
  EXPECT_FALSE(ingress.offer(1, another, another->size()));
}

} // namespace
