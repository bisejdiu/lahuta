/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::array<std::string_view, 4> parts{"besian", "", "sejdiu", "@gmail.com"};
 *   std::vector<std::string_view> valid, empty;
 *   std::partition_copy(parts.begin(), parts.end(), std::back_inserter(valid), std::back_inserter(empty),
 *     [](std::string_view s) { return !s.empty(); });
 *   std::string s; for (auto p : valid) s += p; return s;
 * }();
 *
 */

#ifndef LAHUTA_POOL_FACTORY_HPP
#define LAHUTA_POOL_FACTORY_HPP

#include <atomic>
#include <cstddef>
#include <cstdlib>
#include <memory>
#include <mutex>
#include <vector>

#include "models/pools.hpp"

// clang-format off
namespace lahuta {

template <typename PoolType>
class PoolFactory {
public:

  // Initialize the pool factory with the provided number of threads.
  //
  // Grow-only, idempotent. Never shrinks or clears existing pools.
  // This design is chosen for several reasons:
  // - It's much easier to reason about correctness and safety.
  // - Peak memory usage is equal to the maximum number of items for any file, per thread.
  // - In Python, the same thread may read a single file or use the Pipeline wrapper API
  //   to process multiple files. Object lifetime in Python is not deterministic, and the
  //   grow-only design sidesteps these issues. This is not theoretical. If we clear the
  //   pools below, instead of keeping the capacity, we'll regress.
  //   The test file `test_system_and_topology_from_model_input.py` demonstrates this
  //   behavior (and now guards against it).              - Besian, September 2025
  //
  // Policy: if already >= num_threads, do nothing, otherwise append new pools.
  static void initialize(size_t num_threads) {
    std::lock_guard<std::mutex> lock(grow_mtx_);

    // already initialized to sufficient capacity. Do not shrink or clear (i.e., no pool_.clear(), see above).
    if (pools_.size() >= num_threads) return;

    const size_t old = pools_.size();
    pools_.reserve(num_threads);
    for (size_t i = old; i < num_threads; ++i) {
      pools_.push_back(std::make_unique<PoolType>());
      free_indices_.push_back(i);
    }
  }

  /// Get the pool associated with the current thread.
  static PoolType *get_pool_for_current_thread() {
    struct ThreadLocalPool {
      ThreadLocalPool() {
        std::lock_guard<std::mutex> lock(grow_mtx_);
        if (!free_indices_.empty()) {
          index_ = free_indices_.back();
          free_indices_.pop_back();
        } else {
          index_ = pools_.size();
          pools_.push_back(std::make_unique<PoolType>());
        }
        pool_ = pools_[index_].get();
      }

      ~ThreadLocalPool() {
        if (!factory_alive_.load(std::memory_order_relaxed)) return;

        std::lock_guard<std::mutex> lock(grow_mtx_);
        pools_[index_]->reset();
        free_indices_.push_back(index_);
      }

      PoolType *pool_{};
      size_t index_{};
    };

    thread_local ThreadLocalPool local_pool;
    return local_pool.pool_;
  }

  /// Get a fresh pool for the current thread, resetting its state.
  static PoolType *get_fresh_pool_for_current_thread() {
    PoolType *pool = get_pool_for_current_thread();
    pool->reset();
    return pool;
  }

private:
  inline static std::vector<std::unique_ptr<PoolType>> pools_;
  inline static std::mutex grow_mtx_;
  inline static std::vector<size_t> free_indices_;
  inline static std::atomic<bool> factory_alive_{true};

  struct ShutdownHook {
    ShutdownHook() { std::atexit(&PoolFactory::mark_shutdown); }
  };
  inline static ShutdownHook shutdown_hook_{};

  static void mark_shutdown() { factory_alive_.store(false, std::memory_order_relaxed); }
};

// Using aliases is safer, because referencing the template directly may lead to multiple independent instantiations
using AtomPoolFactory = PoolFactory<AtomPool>;
using BondPoolFactory = PoolFactory<BondPool>;
using InfoPoolFactory = PoolFactory<InfoPool>;

} // namespace lahuta

#endif // LAHUTA_POOL_FACTORY_HPP
