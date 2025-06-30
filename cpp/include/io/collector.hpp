#ifndef LAHUTA_IO_COLLECTOR_HPP
#define LAHUTA_IO_COLLECTOR_HPP

#include <atomic>
#include <condition_variable>
#include <cstddef>
#include <mutex>

#include "pipeline/core/emitter.hpp"

//
// A global counter when we have multi collectors is an anti-pattern:
//  - one mistake in one collector can block all others (this actually makes it good for forcing correctness)
//  - miscounting size may lead to underflow/overflow
//  - no way to know which collector is blocked
// Instead, each collector could have its own counter, or pass the counter via DI in the form of a controller
// object (would make it easier to test).  - Besian, May 16, 2025
//

// clang-format off
namespace lahuta {

// global backpressure counters (see fix note above)
inline std::atomic<std::size_t> g_bytes_in_flight{0};
inline std::size_t HIGH_WATER = 512ul * 1024 * 1024; // 512MB
inline std::size_t LOW_WATER  = 300ul * 1024 * 1024; // 300MB // NOTE: not used

template<typename T, typename Policy>
class Collector final : public pipeline::IEmitter<T> {
public:
  using ptr_type = typename pipeline::IEmitter<T>::ptr_type;
  using input_type = ptr_type;

  template<typename... Args>
  explicit Collector(Args&&... args)
    : policy_(std::forward<Args>(args)...) {}

  //
  // Both versions are inefficient. At some point we may want to revisit these, but I suspect
  // larger refactor to upstream and downstream code is needed to make this work.
  //

  // void emit(const std::shared_ptr<const T>& item) {
  //   std::unique_lock lk(m_);
  //   wait_for_space(lk);

  //   const auto delta = policy_.append_size(*item);
  //   policy_.append(*item);
  //   g_bytes_in_flight += delta;

  //   if (policy_.needs_flush()) {
  //     auto freed = policy_.flush();
  //     g_bytes_in_flight -= freed;
  //     cv_.notify_all();
  //   }
  // }

  void emit(ptr_type item) override {
    const std::size_t delta = policy_.append_size(item);

    std::unique_lock lk(m_);
    for (;;) {
      if (g_bytes_in_flight + delta <= HIGH_WATER) break;

      const std::size_t freed = policy_.flush(); // flushing to a target value, LOW_WATER, could be better
      g_bytes_in_flight -= freed;
    }

    policy_.append(std::move(item));
    g_bytes_in_flight += delta;

    if (policy_.needs_flush()) {
      const auto freed = policy_.flush();
      g_bytes_in_flight -= freed;
    }
    lk.unlock();
    cv_.notify_all();
  }

  void finish() {
    std::scoped_lock lk(m_);
    const auto freed = policy_.finish();
    g_bytes_in_flight -= freed;
    cv_.notify_all();
  }

  decltype(auto) result() const { return policy_.result(); }

private:
  void wait_for_space(std::unique_lock<std::mutex>& lk) {
    cv_.wait(lk, []{return g_bytes_in_flight.load(std::memory_order_relaxed) < HIGH_WATER;});
  }

  Policy                policy_;
  std::mutex            m_;
  std::condition_variable cv_;
};
} // namespace lahuta

#endif // LAHUTA_IO_COLLECTOR_HPP
