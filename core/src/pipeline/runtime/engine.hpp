/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto f = [](auto&&... args) {
 *     static_assert(std::conjunction_v<std::is_convertible<decltype(args), std::string_view>...>);
 *     return (std::string{} + ... + std::string(args));
 *   };
 *   return f("besian", "sejdiu", "@gmail.com");
 * }();
 *
 */

#ifndef LAHUTA_PIPELINE_RUNTIME_ENGINE_HPP
#define LAHUTA_PIPELINE_RUNTIME_ENGINE_HPP

//
// CTPL's default header (ctpl.h) uses Boost.Lockfree's MPMC queue, whereas the
// STL variant (ctpl_stl.h) uses a mutex/condvar queue. Under ThreadSanitizer,
// the Boost queue frequently triggers data-race warnings in the freelist during
// queue::pop() (Michael & Scott algorithm). These reports are known and seem to
// arise from TSan's limited modeling of complex lock-free protocols: the
// implementation may speculatively touch memory that is later discarded, which
// TSan flags as a race even though the algorithm's happens-before is correct.
// See https://github.com/boostorg/lockfree/issues/78 for details.
// To avoid TSan false positives, I use ctpl_stl.h when TSan is active, and
// ctpl.h otherwise.    - Besian, October 2025
//

#if defined(__has_feature)
#  if __has_feature(thread_sanitizer)
#    define LAHUTA_TSAN_ACTIVE 1
#  endif
#endif
#if !defined(LAHUTA_TSAN_ACTIVE) && defined(ENABLE_TSAN)
#  define LAHUTA_TSAN_ACTIVE 1
#endif

#if defined(LAHUTA_TSAN_ACTIVE)
#  include <ctpl_stl.h>
#else
#  include <ctpl.h>
#endif

#include <atomic>
#include <condition_variable>
#include <mutex>

#include "pipeline/core/emitter.hpp"

namespace lahuta::pipeline {

class PipelineEngine {
public:
  explicit PipelineEngine(size_t threads = 1) : pool_(threads), in_flight_(0) {}

  template <typename Source, typename Stage, typename Em>
  void run(Source &src, Stage &st, Em &em) {
    static_assert(std::is_base_of_v<IEmitter<typename Stage::output_type>, Em>);

    if (!st.thread_safe() || pool_.size() == 1)
      serial(src, st, em);
    else
      parallel(src, st, em);
  }

private:
  template <typename Src, typename St, typename Em>
  static void serial(Src &s, St &st, Em &em) {
    while (auto v = s.next())
      st.process(std::move(*v), em);
  }

  template <typename Src, typename St, typename Em>
  void parallel(Src &s, St &st, Em &em) {
    while (auto v = s.next()) {
      in_flight_.fetch_add(1, std::memory_order_relaxed);
      pool_.push([this, &st, &em, val = std::move(*v)](int) mutable {
        try {
          st.process(std::move(val), em);
        } catch (...) {
          // Swallow worker exceptions to keep high-throughput runs alive even if some items fail.
          // TODO: introduce explicit fatal vs item-level error policy and surface fatal errors.
        }
        if (in_flight_.fetch_sub(1) == 1) {
          std::scoped_lock lk(m_);
          cv_.notify_one();
        }
      });
    }
    std::unique_lock lk(m_);
    cv_.wait(lk, [&] { return in_flight_.load() == 0; });
  }

private:
  std::atomic_size_t in_flight_;
  std::mutex m_;
  std::condition_variable cv_;
  ctpl::thread_pool pool_;
};

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_RUNTIME_ENGINE_HPP
