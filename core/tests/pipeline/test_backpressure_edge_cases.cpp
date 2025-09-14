// Additional tests for dynamic backpressure covering:
// - emit behavior after multiplexer is closed
// - rotation by record count
// - non-required sink failure behavior during close_and_flush
#include <gtest/gtest.h>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "io/sinks/memory.hpp"
#include "io/sinks/sharded_ndjson.hpp"
#include "pipeline/dynamic/sink_iface.hpp"
#include <pipeline/dynamic/channel_multiplexer.hpp>

using namespace lahuta::pipeline::dynamic;

// clang-format off
namespace {
class ThrowingSink : public IDynamicSink {
public:
  void write(EmissionView) override { throw std::runtime_error("boom"); }
};
} // namespace

// Emit after close should be ignored
// After a successful close_and_flush, further emit() calls must be no-ops and
// must not append to sink outputs. Second close is also a no-op.
TEST(DynamicBackpressure_EdgeCases, EmitIgnoredAfterClose) {
  ChannelMultiplexer mux;
  auto mem = std::make_shared<MemorySink>();
  mux.connect("x", mem);
  mux.emit(Emission{"x", "one"});
  ASSERT_NO_THROW(mux.close_and_flush(std::chrono::seconds(1)));
  auto v1 = mem->result();
  ASSERT_EQ(v1.size(), 1u);

  // Post-close emit should be ignored
  mux.emit(Emission{"x", "two"});
  // Second close is a no-op when already closed
  ASSERT_NO_THROW(mux.close_and_flush(std::chrono::milliseconds(10)));
  auto v2 = mem->result();
  EXPECT_EQ(v2.size(), 1u);
}

// Count-based rotation produces expected number of shards and line counts
// With shard_size=2 and 5 records, we expect shards with 2,2,1 lines.
TEST(DynamicBackpressure_EdgeCases, ShardedSinkRotatesByCount) {
  namespace fs = std::filesystem;
  fs::path dir = fs::temp_directory_path() / fs::path("lahuta_test_shards_count");
  std::error_code ec; fs::create_directories(dir, ec);
  for (auto &p : fs::directory_iterator(dir, ec)) { (void)fs::remove(p.path(), ec); }

  auto sink = std::make_shared<ShardedNdjsonSink>(dir.string(), /*shard_size=*/2);
  ChannelMultiplexer mux;
  mux.connect("s", sink);

  mux.emit(Emission{"s", "a"});
  mux.emit(Emission{"s", "b"});
  mux.emit(Emission{"s", "c"});
  mux.emit(Emission{"s", "d"});
  mux.emit(Emission{"s", "e"});

  ASSERT_NO_THROW(mux.close_and_flush(std::chrono::seconds(2)));

  auto files = sink->files();
  ASSERT_EQ(files.size(), 3u);

  // Expect 2,2,1 lines across shards
  std::vector<size_t> counts;
  for (const auto& f : files) {
    std::ifstream in(f);
    ASSERT_TRUE(in.good());
    size_t n = 0; std::string line;
    while (std::getline(in, line)) ++n;
    counts.push_back(n);
  }
  ASSERT_EQ(counts.size(), 3u);
  EXPECT_EQ(counts[0], 2u);
  EXPECT_EQ(counts[1], 2u);
  EXPECT_EQ(counts[2], 1u);

  for (const auto& f : files) { (void)fs::remove(f, ec); }
  (void)fs::remove(dir, ec);
}

// Writer throws but required=false: close_and_flush should succeed
// A non-required sink that fails in write() must not abort the pipeline. The
// multiplexer should still report success at shutdown.
TEST(DynamicBackpressure_EdgeCases, NonRequiredSinkFailureDoesNotAbort) {
  ChannelMultiplexer mux;
  auto bad = std::make_shared<ThrowingSink>();
  BackpressureConfig cfg; cfg.required = false;
  mux.connect("err", bad, cfg);
  mux.emit(Emission{"err", "x"});
  EXPECT_NO_THROW(mux.close_and_flush(std::chrono::seconds(2)));
}
