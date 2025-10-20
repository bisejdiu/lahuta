// Exercises the dynamic pipeline's bounded backpressure implementation.
//
// These focus on the core invariants:
// - bounded queues, by messages and bytes
// - configurable full-buffer policies (block/drop/coalesce)
// - per-sink isolation (slow sink does not stall others)
// - deterministic shutdown semantics (drain/timeout/failure propagation)
// - file rotation policies (by bytes)
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <filesystem>
#include <fstream>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include <gtest/gtest.h>
#include <pipeline/dynamic/backpressure.hpp>
#include <pipeline/dynamic/channel_multiplexer.hpp>

#include "io/sinks/memory.hpp"
#include "io/sinks/sharded_ndjson.hpp"
#include "pipeline/dynamic/sink_iface.hpp"

using namespace lahuta::pipeline::dynamic;

// clang-format off
namespace {

// Build a QueueNode from text and channel id
static QueueNode make_node(uint32_t ch_id, const std::string& s) {
  auto buf = std::make_shared<std::string>(s);
  return QueueNode{ch_id, buf, buf->size()};
}

// Slow writer that sleeps N milliseconds per write
class SleepySink : public IDynamicSink {
public:
  explicit SleepySink(int sleep_ms) : sleep_ms_(sleep_ms) {}

  void write(EmissionView e) override {
    ++writes_;
    if (sleep_ms_ > 0) std::this_thread::sleep_for(std::chrono::milliseconds(sleep_ms_));
    last_size_.store(e.payload.size(), std::memory_order_relaxed);
  }
  void flush() override {}
  void close() override {}

  uint64_t  writes() const { return writes_   .load(std::memory_order_relaxed); }
  size_t last_size() const { return last_size_.load(std::memory_order_relaxed); }

private:
  int sleep_ms_;
  std::atomic<uint64_t> writes_{0};
  std::atomic<size_t> last_size_{0};
};

// Blocks indefinitely in write() until released
class BlockingSink : public IDynamicSink {
public:
  void write(EmissionView) override {
    std::unique_lock<std::mutex> lk(m_);
    cv_.wait(lk, [&]{ return released_; });
  }
  void flush() override {}
  void close() override {}

  void release() {
    std::lock_guard<std::mutex> lk(m_);
    released_ = true;
    cv_.notify_all();
  }

private:
  std::mutex m_;
  std::condition_variable cv_;
  bool released_ = false;
};

// Throwing sink that raises in write()
class ThrowingSink : public IDynamicSink {
public:
  void write(EmissionView) override { throw std::runtime_error("write failed"); }
  void flush() override {}
  void close() override {}
};

} // namespace

// Blocking under full. Unit-test the BoundedQueue.
// Verifies that a producer attempting to enqueue to a full queue blocks until
// the consumer drains space, accumulates measurable stall time, and ultimately
// succeeds without drops under OnFull::Block.
TEST(DynamicBackpressure, BoundedQueueBlocksWhenFullUntilDrained) {
  BoundedQueue q(/*max_msgs=*/1, /*max_bytes=*/1024 * 1024);
  BackpressureConfig cfg;
  cfg.on_full = OnFull::Block;
  cfg.offer_wait_slice = std::chrono::milliseconds(10);

  std::atomic<uint64_t> stall_ns{0};
  std::atomic<uint64_t> drops{0};

  // Fill one slot
  ASSERT_TRUE(q.offer(make_node(1, "A"), cfg, &stall_ns, &drops));

  // Start a thread that will attempt to offer and must block until consumer drains
  std::atomic<bool> offered{false};
  auto producer = std::thread([&]{
    bool ok = q.offer(make_node(1, "B"), cfg, &stall_ns, &drops);
    offered.store(ok, std::memory_order_relaxed);
  });

  // Sleep to ensure the producer experiences at least one timed wait
  std::this_thread::sleep_for(std::chrono::milliseconds(50));

  // Drain queue: pop_batch should release capacity
  std::vector<QueueNode> batch;
  ASSERT_TRUE(q.pop_batch(batch, /*max_msgs=*/8, /*max_bytes=*/1024 * 1024));
  ASSERT_EQ(batch.size(), 1u);

  producer.join();
  ASSERT_TRUE(offered.load());
  // We should have accumulated at least ~ one wait period worth of stall time
  EXPECT_GE(stall_ns.load(), static_cast<uint64_t>(5'000'000)); // >= 5ms
  EXPECT_EQ(drops.load(), 0u);
}

// DropLatest / DropOldest policies
// Ensures that when full:
// - DropLatest rejects the new item and increments drop counters
// - DropOldest evicts the oldest item to make room so the newest is preserved
TEST(DynamicBackpressure, DropPoliciesBehaveAsExpected) {
  // DropLatest
  {
    BoundedQueue q(1, 1024 * 1024);
    BackpressureConfig cfg; cfg.on_full = OnFull::DropLatest;
    std::atomic<uint64_t> stall{0}, drops{0};
    ASSERT_TRUE(q.offer(make_node(1, "first"), cfg, &stall, &drops));
    EXPECT_FALSE(q.offer(make_node(1, "second"), cfg, &stall, &drops));
    EXPECT_EQ(drops.load(), 1u);
  }

  // DropOldest ensures newest is kept
  {
    BoundedQueue q(1, 1024 * 1024);
    BackpressureConfig cfg; cfg.on_full = OnFull::DropOldest;
    std::atomic<uint64_t> stall{0}, drops{0};
    ASSERT_TRUE(q.offer(make_node(1, "old"), cfg, &stall, &drops));
    EXPECT_TRUE(q.offer(make_node(1, "new"), cfg, &stall, &drops));
    std::vector<QueueNode> batch;
    ASSERT_TRUE(q.pop_batch(batch, /*max_msgs=*/8, /*max_bytes=*/1024 * 1024));
    ASSERT_EQ(batch.size(), 1u);
    EXPECT_EQ(batch.front().size, std::string("new").size());
    EXPECT_EQ(drops.load(), 1u);
  }
}

TEST(DynamicBackpressure, DefaultConfigSetterAppliesToNewConnections) {
  auto original = get_default_backpressure_config();
  struct Reset {
    BackpressureConfig cfg;
    ~Reset() { set_default_backpressure_config(cfg); }
  } guard{original};

  BackpressureConfig custom = original;
  custom.max_queue_bytes = 512;
  custom.max_batch_bytes = 512;
  custom.on_full = OnFull::DropLatest;
  custom.required = false;
  set_default_backpressure_config(custom);

  ChannelMultiplexer mux;
  auto sink = std::make_shared<MemorySink>();
  mux.connect("test", sink);

  Emission e;
  e.channel = "test";
  e.payload.assign(2048, 'x');
  mux.emit(std::move(e));

  mux.close_and_flush(std::chrono::milliseconds(100));

  auto stats = mux.stats();
  ASSERT_EQ(stats.size(), 1u);
  EXPECT_EQ(stats[0].drops, 1u);
  EXPECT_TRUE(sink->result().empty());
}

// Fairness (isolation): slow sink must not block fast sink, and fast sink should receive all items.
// Connect a slow sink (sleeping writer) and a fast sink to the same channel.
// Emit many records and confirm the fast sink receives all while the slow sink
// may lag/drop according to its policy. This validates per-sink isolation.
TEST(DynamicBackpressure, SlowSinkDoesNotBlockFastSink) {
  ChannelMultiplexer mux;
  auto slow = std::make_shared<SleepySink>(1); // 1ms per write
  auto fast = std::make_shared<SleepySink>(0);

  BackpressureConfig slow_cfg; slow_cfg.on_full = OnFull::DropLatest; slow_cfg.max_queue_msgs = 64;
  BackpressureConfig fast_cfg; fast_cfg.on_full = OnFull::Block;      fast_cfg.max_queue_msgs = 64;

  mux.connect("ch", slow, slow_cfg);
  mux.connect("ch", fast, fast_cfg);

  const int N = 5000;
  for (int i = 0; i < N; ++i) {
    mux.emit(Emission{"ch", std::to_string(i)});
  }
  ASSERT_NO_THROW(mux.close_and_flush(std::chrono::seconds(5)));

  // Fast sink must see all N, slow sink may drop
  EXPECT_EQ(fast->writes(), static_cast<uint64_t>(N));
  EXPECT_LE(slow->writes(), static_cast<uint64_t>(N));
}

// Memory ceiling
// Enforces the per-sink max_queue_bytes budget by rejecting enqueues that
// would exceed the byte cap (DropLatest policy), and increments drop counts.
TEST(DynamicBackpressure, MaxQueueBytesEnforced) {
  BoundedQueue q(/*max_msgs=*/10, /*max_bytes=*/100);
  BackpressureConfig cfg; cfg.on_full = OnFull::DropLatest;
  std::atomic<uint64_t> stall{0}, drops{0};

  std::string big(80, 'X');
  ASSERT_TRUE(q.offer(make_node(1, big), cfg, &stall, &drops));
  // Another big entry would exceed 100 bytes
  EXPECT_FALSE(q.offer(make_node(1, big), cfg, &stall, &drops));
  EXPECT_EQ(drops.load(), 1u);
}

// Oversize relative to batch budget should still drain in a single pop_batch() iteration.
TEST(DynamicBackpressure, PopBatchMakesProgressForOversizedEntry) {
  BoundedQueue q(/*max_msgs=*/8, /*max_bytes=*/1024);
  BackpressureConfig cfg;
  cfg.max_batch_msgs = 4;
  cfg.max_batch_bytes = 64;
  std::atomic<uint64_t> stall{0}, drops{0};

  std::string big(128, 'Z');
  ASSERT_TRUE(q.offer(make_node(1, big), cfg, &stall, &drops));

  std::vector<QueueNode> batch;
  ASSERT_TRUE(q.pop_batch(batch, cfg.max_batch_msgs, cfg.max_batch_bytes));
  ASSERT_EQ(batch.size(), 1u);
  EXPECT_EQ(batch.front().size, big.size());
}

// Items larger than the queue byte budget must be rejected immediately to avoid blocking forever.
TEST(DynamicBackpressure, OversizedItemRejectedWhenExceedingQueueBytes) {
  BoundedQueue q(/*max_msgs=*/4, /*max_bytes=*/32);
  BackpressureConfig cfg; cfg.on_full = OnFull::Block;
  std::atomic<uint64_t> stall{0}, drops{0};

  std::string huge(64, 'H');
  EXPECT_FALSE(q.offer(make_node(1, huge), cfg, &stall, &drops));
  EXPECT_EQ(drops.load(), 1u);
  EXPECT_EQ(stall.load(), 0u);
  EXPECT_EQ(q.size_msgs(), 0u);
}

// Clean shutdown drains and matches counts
// Emits N records and confirms close_and_flush drains the queue and writers,
// with the sink observing exactly N outputs (no losses beyond policy).
TEST(DynamicBackpressure, CleanShutdownDrains) {
  ChannelMultiplexer mux;
  auto mem = std::make_shared<MemorySink>();
  mux.connect("contacts", mem);

  const int N = 1000;
  for (int i = 0; i < N; ++i) {
    mux.emit(Emission{"contacts", std::string("rec_") + std::to_string(i)});
  }
  ASSERT_NO_THROW(mux.close_and_flush(std::chrono::seconds(5)));
  auto out = mem->result();
  EXPECT_EQ(static_cast<int>(out.size()), N);
}

// Timeout on shutdown returns failure for required sink
// A sink whose write() blocks causes close_and_flush to time out for required
// sinks. The test asserts an error is surfaced, and uses a watchdog to release
// the sink to avoid dangling threads.
TEST(DynamicBackpressure, TimeoutOnShutdownFailsForRequiredSink) {
  ChannelMultiplexer mux;
  auto blk = std::make_shared<BlockingSink>();
  BackpressureConfig cfg; cfg.required = true; cfg.max_queue_msgs = 1;
  mux.connect("a", blk, cfg);

  // Emit one record so writer blocks in write()
  mux.emit(Emission{"a", "x"});

  // Watchdog releases the sink after a delay to ensure the writer eventually exits
  std::thread watchdog([blk]{
    std::this_thread::sleep_for(std::chrono::milliseconds(200));
    blk->release();
  });

  bool threw = false;
  try {
    mux.close_and_flush(std::chrono::milliseconds(50));
  } catch (...) {
    threw = true;
  }
  EXPECT_TRUE(threw);
  watchdog.join();
}

// Rotation by bytes: verify shard sizes respect threshold (+ last record allowance)
// Configures a bytes threshold and emits fixed-size records, then verifies that
// shard files are produced and each file size respects the limit, allowing for
// the final record boundary.
TEST(DynamicBackpressure, ShardedSinkRotatesByBytes) {
  namespace fs = std::filesystem;
  fs::path tmp = fs::temp_directory_path() / fs::path("lahuta_test_shards");
  std::error_code ec; fs::create_directories(tmp, ec);
  // Clean up any previous files
  for (auto &p : fs::directory_iterator(tmp, ec)) { (void)fs::remove(p.path(), ec); }

  auto sink = std::make_shared<ShardedNdjsonSink>(tmp.string(), /*shard_size=*/1000000, /*max_shard_bytes=*/1024);
  ChannelMultiplexer mux;
  mux.connect("s", sink);

  // Emit records of 200 bytes to force rotations
  const std::string rec(200, 'A');
  for (int i = 0; i < 20; ++i) mux.emit(Emission{"s", rec});
  ASSERT_NO_THROW(mux.close_and_flush(std::chrono::seconds(5)));

  auto files = sink->files();
  ASSERT_GT(files.size(), 1u);

  // Each file size should be <= threshold + one record (allowance at boundary)
  for (const auto& f : files) {
    std::ifstream in(f, std::ios::binary | std::ios::ate);
    ASSERT_TRUE(in.good());
    auto sz = static_cast<size_t>(in.tellg());
    EXPECT_LE(sz, static_cast<size_t>(1024 + rec.size() + 1));
  }

  // Cleanup
  for (const auto& f : files) { (void)fs::remove(f, ec); }
  (void)fs::remove(tmp, ec);
}

// Error propagation from sink writer
// A sink that throws inside write() should cause close_and_flush to fail when
// the sink is required=true, validating failure propagation.
TEST(DynamicBackpressure, WriterExceptionPropagatesOnCloseForRequired) {
  ChannelMultiplexer mux;
  auto bad = std::make_shared<ThrowingSink>();
  BackpressureConfig cfg; cfg.required = true;
  mux.connect("err", bad, cfg);
  mux.emit(Emission{"err", "x"});
  EXPECT_ANY_THROW(mux.close_and_flush(std::chrono::seconds(2)));
}
