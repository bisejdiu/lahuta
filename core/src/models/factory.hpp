#ifndef LAHUTA_POOL_FACTORY_HPP
#define LAHUTA_POOL_FACTORY_HPP

#include <atomic>
#include <cassert>
#include <cstddef>
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
      if (pools_.size() >= num_threads) {
          return; // already initialized to sufficient capacity. Do not shrink or clear (i.e., no pool_.clear(), see above).
      }
      const size_t old = pools_.size();
      pools_.reserve(num_threads);
      for (size_t i = old; i < num_threads; ++i) {
          pools_.push_back(std::make_unique<PoolType>());
      }
  }

  /// Get a pool for a specific thread ID.
  static PoolType* get_pool_by_thread_id(size_t thread_id) {
      assert(thread_id < pools_.size() && "Thread ID out of range");
      return pools_[thread_id].get();
  }

  /// Get the pool associated with the current thread.
  static PoolType* get_pool_for_current_thread() {
      static std::atomic<size_t> next_id{0};
      static thread_local size_t thread_index = next_id.fetch_add(1, std::memory_order_relaxed);

      // Self-initialize on first touch for this thread.
      if (thread_index >= pools_.size()) {
          std::lock_guard<std::mutex> lock(grow_mtx_);
          if (thread_index >= pools_.size()) {
              const size_t old = pools_.size();
              pools_.reserve(thread_index + 1);
              for (size_t i = old; i <= thread_index; ++i) {
                  pools_.push_back(std::make_unique<PoolType>());
              }
          }
      }
      return pools_[thread_index].get();
  }

  /// Get a fresh pool for the current thread, resetting its state.
  static PoolType* get_fresh_pool_for_current_thread() {
      PoolType* pool = get_pool_for_current_thread();
      pool->reset();
      return pool;
  }

private:
  inline static std::vector<std::unique_ptr<PoolType>> pools_;
  inline static std::mutex grow_mtx_;
};

// Using aliases is safer, because referencing the template directly could lead to
// multiple independent instantiations
using AtomPoolFactory = PoolFactory<AtomPool>;
using BondPoolFactory = PoolFactory<BondPool>;
using InfoPoolFactory = PoolFactory<InfoPool>;

} // namespace lahuta

#endif // LAHUTA_POOL_FACTORY_HPP
