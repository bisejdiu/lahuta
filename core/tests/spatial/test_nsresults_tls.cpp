#include <atomic>
#include <cstdint>
#include <thread>
#include <vector>

#include <gtest/gtest.h>

#include "spatial/nsresults_tls.hpp"

// clang-format off
namespace lahuta {

TEST(NSResultsTLS, ReusesSameInstanceWithinThread) {
  NSResults* first_ptr = nullptr;
  std::size_t first_capacity = 0;

  {
    TlsResultsScope scope;
    auto& results = scope.results();
    results.reserve(128);
    results.add_neighbors(1, 2, 3.0f);

    EXPECT_EQ(results.size(), 1);
    first_ptr = &results;
    first_capacity = results.get_pairs().capacity();
    EXPECT_GE(first_capacity, 1u);
  }

  {
    TlsResultsScope scope;
    auto& results = scope.results();
    EXPECT_EQ(&results, first_ptr);
    EXPECT_EQ(results.size(), 0u);
    EXPECT_GE(results.get_pairs().capacity(), first_capacity);
  }
}

TEST(NSResultsTLS, ProvidesDistinctInstancesPerThread) {
  NSResults* main_thread_ptr = nullptr;
  {
    TlsResultsScope scope;
    main_thread_ptr = &scope.results();
  }

  std::atomic<NSResults*> worker_ptr{nullptr};
  std::thread worker([&]() {
    TlsResultsScope scope;
    worker_ptr.store(&scope.results(), std::memory_order_relaxed);
    EXPECT_EQ(scope.results().size(), 0u);
  });
  worker.join();

  ASSERT_NE(worker_ptr.load(std::memory_order_relaxed), nullptr);
  EXPECT_NE(worker_ptr.load(std::memory_order_relaxed), main_thread_ptr);
}

TEST(NSResultsTLS, CapacityResetWhenExceedingSoftMax) {
  NSResults* first_instance = nullptr;
  std::size_t original_capacity = 0;
  std::size_t original_size = 0;

  // First, create an instance and grow it beyond the soft max
  {
    TlsResultsScope scope;
    auto& results = scope.results();

    results.reserve(NSRESULTS_TLS_SOFT_MAX + 1000);
    results.add_neighbors(1, 2, 3.0f);

    first_instance = &results;
    original_capacity = results.get_pairs().capacity();
    original_size = results.size();
    EXPECT_GT(original_capacity, NSRESULTS_TLS_SOFT_MAX);
    EXPECT_EQ(original_size, 1u);
  }

  // Second scope should get the same instance but with reset contents
  {
    TlsResultsScope scope;
    auto& results = scope.results();

    EXPECT_EQ(&results, first_instance); // Should be the same instance address
    EXPECT_EQ(results.size(), 0u);       // But contents should be reset
    results.reserve(10);                 // Capacity should be reset to minimal
    EXPECT_LT(results.get_pairs().capacity(), original_capacity / 2);
  }
}

TEST(NSResultsTLS, PreservesCapacityWhenBelowSoftMax) {
  NSResults* first_instance = nullptr;
  std::size_t small_capacity = 0;

  {
    TlsResultsScope scope;
    auto& results = scope.results();
    results.reserve(1000);  // below 2M soft max
    results.add_neighbors(1, 2, 3.0f);

    first_instance = &results;
    small_capacity = results.get_pairs().capacity();
    EXPECT_GT(small_capacity, 0u);
    EXPECT_LT(small_capacity, NSRESULTS_TLS_SOFT_MAX);
  }

  {
    TlsResultsScope scope;
    auto& results = scope.results();

    EXPECT_EQ(&results, first_instance); // Should be the same instance
    EXPECT_EQ(results.size(), 0u);       // Should start empty but preserve capacity
    EXPECT_GE(results.get_pairs().capacity(), small_capacity);
  }
}

TEST(NSResultsTLS, MultipleConsecutiveLargeAllocations) {
  NSResults* first_instance = nullptr;
  std::size_t large_capacity = 0;

  {
    TlsResultsScope scope;
    auto& results = scope.results();

    results.reserve(NSRESULTS_TLS_SOFT_MAX + 10000); // Trigger capacity reset

    for (int j = 0; j < 100; ++j) {
      results.add_neighbors(j, j + 1, static_cast<float>(j));
    }

    first_instance = &results;
    large_capacity = results.get_pairs().capacity();
    EXPECT_GT(large_capacity, NSRESULTS_TLS_SOFT_MAX);
  }

  {
    TlsResultsScope scope;
    auto& results = scope.results();

    EXPECT_EQ(&results, first_instance); // Same instance address (thread-local reuse)
    EXPECT_EQ(results.size(), 0u);       // But contents should be cleared
    results.reserve(100);                // And capacity should be reset to minimal
    EXPECT_LT(results.get_pairs().capacity(), large_capacity / 2);
  }

  {
    TlsResultsScope scope;
    auto& results = scope.results();

    EXPECT_EQ(&results, first_instance);
    EXPECT_EQ(results.size(), 0u);
  }
}

TEST(NSResultsTLS, MultipleWorkersIsolationAndSoftMaxReset) {
  constexpr int kNumThreads = 8;

  struct ThreadStats {
    std::uintptr_t addr_first = 0;
    std::uintptr_t addr_second = 0;
    std::size_t cap1_pairs = 0;
    std::size_t cap1_dists = 0;
    std::size_t cap2_pairs = 0;
    std::size_t cap2_dists = 0;
    std::size_t size1 = 0;
    std::size_t size2 = 0;
  };

  std::vector<std::thread> workers;
  std::vector<ThreadStats> stats(kNumThreads);
  std::atomic<int> ready{0};
  std::atomic<bool> go{false};

  workers.reserve(kNumThreads);
  for (int i = 0; i < kNumThreads; ++i) {
    workers.emplace_back([i, &stats, &ready, &go]() {
      ThreadStats local;

      {
        TlsResultsScope scope;
        auto &results = scope.results();

        results.reserve(NSRESULTS_TLS_SOFT_MAX + 5000);
        for (int j = 0; j < 256; ++j) {
          results.add_neighbors(i * 1000 + j, i * 1000 + j + 1, static_cast<float>(j));
        }

        local.addr_first = reinterpret_cast<std::uintptr_t>(&results);
        local.cap1_pairs = results.get_pairs().capacity();
        local.cap1_dists = results.get_distances().capacity();
        local.size1 = results.size();

        EXPECT_GT(std::max(local.cap1_pairs, local.cap1_dists), NSRESULTS_TLS_SOFT_MAX);
        EXPECT_EQ(local.size1, 256u);

        // Simple barrier so all threads overlap while holding their TLS instance
        int prev = ready.fetch_add(1, std::memory_order_acq_rel) + 1;
        if (prev == kNumThreads) {
          go.store(true, std::memory_order_release);
        } else {
          while (!go.load(std::memory_order_acquire)) {
          }
        }

        EXPECT_EQ(results.size(), 256u); // No cross-thread interference changed our size
      }

      {
        TlsResultsScope scope;
        auto &results = scope.results();

        local.addr_second = reinterpret_cast<std::uintptr_t>(&results);
        local.size2 = results.size();

        // Instance should have been reset
        EXPECT_EQ(local.addr_second, local.addr_first);
        EXPECT_EQ(local.size2, 0u);

        results.reserve(128);
        local.cap2_pairs = results.get_pairs().capacity();
        local.cap2_dists = results.get_distances().capacity();

        const std::size_t cap1_max = std::max(local.cap1_pairs, local.cap1_dists);
        const std::size_t cap2_max = std::max(local.cap2_pairs, local.cap2_dists);
        EXPECT_LE(cap2_max, NSRESULTS_TLS_SOFT_MAX);
        EXPECT_LT(cap2_max, cap1_max);
      }

      stats[i] = local;
    });
  }

  for (auto &t : workers) {
    t.join();
  }

  // Verify reuse within each thread
  for (const auto &s : stats) {
    ASSERT_NE(s.addr_first, 0u);
    ASSERT_NE(s.addr_second, 0u);
    EXPECT_EQ(s.addr_first, s.addr_second);
  }
}

} // namespace lahuta
