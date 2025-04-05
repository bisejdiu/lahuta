#ifndef LAHUTA_POOL_FACTORY_HPP
#define LAHUTA_POOL_FACTORY_HPP

#include "models/pools.hpp"
#include <cstddef>
#include <memory>
#include <vector>

namespace lahuta {

template <typename PoolType>
class PoolFactory {
public:
    /// Initialize the pool factory with a specified number of threads.
    static void initialize(size_t num_threads) {
        pools_.clear();
        pools_.reserve(num_threads);
        for (size_t i = 0; i < num_threads; ++i) {
            pools_.push_back(std::make_unique<PoolType>());
        }
    }

    /// Get a pool for a specific thread ID.
    static PoolType* getPoolByThreadId(size_t thread_id) {
        assert(thread_id < pools_.size() && "Thread ID out of range");
        return pools_[thread_id].get();
    }

    /// Get the pool associated with the current thread.
    static PoolType* getPoolForCurrentThread() {
        static std::atomic<size_t> next_id{0};
        static thread_local size_t thread_index = next_id++;

        assert(thread_index < pools_.size() && "Thread index exceeds initialized pool count");
        return getPoolByThreadId(thread_index);
    }

    /// Get a fresh pool for the current thread, resetting its state.
    static PoolType* getFreshPoolForCurrentThread() {
        PoolType* pool = getPoolForCurrentThread();
        pool->reset();
        return pool;
    }

private:
    inline static std::vector<std::unique_ptr<PoolType>> pools_;
};

// Using aliases is safer, because referencing the template directly could lead to
// multiple independent instantiations
using AtomPoolFactory = PoolFactory<AtomPool>;
using BondPoolFactory = PoolFactory<BondPool>;
using InfoPoolFactory = PoolFactory<InfoPool>;

} // namespace lahuta

#endif // LAHUTA_POOL_FACTORY_HPP
