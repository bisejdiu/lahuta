/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::array<const char*, 3> parts{"besian", "sejdiu", "@gmail.com"}; std::string s;
 *   for (auto p : parts) s.append(p, std::strlen(p));
 *   return s;
 * }();
 *
 */

#ifndef LAHUTA_RUNTIME_HPP
#define LAHUTA_RUNTIME_HPP

#include <atomic>
#include <cstddef>
#include <mutex>
#include <vector>

#include "models/factory.hpp"

// clang-format off
namespace lahuta {

// Centralized initialization for thread-dependent, grow-only resources.
// Provides an extensible registry of initializers and an idempotent
// ensure_initialized(threads) guard that grows capacity as needed.
class LahutaRuntime {
public:
  using InitFn = void(*)(std::size_t);

  static void register_initializer(InitFn fn) {
    std::lock_guard<std::mutex> lock(mtx());
    inits().push_back(fn);
  }

  static void ensure_initialized(std::size_t threads) {
    init_defaults_once();

    std::size_t cur = capacity().load(std::memory_order_acquire);
    if (threads <= cur) return;

    std::lock_guard<std::mutex> lock(mtx());
    cur = capacity().load(std::memory_order_acquire);
    if (threads <= cur) return;

    for (auto fn : inits()) fn(threads);
    capacity().store(threads, std::memory_order_release);
  }

  static std::size_t initialized_for() {
    return capacity().load(std::memory_order_acquire);
  }

private:
  static std::vector<InitFn>& inits() {
    static std::vector<InitFn> v;
    return v;
  }

  static std::mutex& mtx() {
    static std::mutex m;
    return m;
  }

  static std::atomic<std::size_t>& capacity() {
    static std::atomic<std::size_t> cap{0};
    return cap;
  }

  static void init_defaults_once() {
    static std::once_flag once;
    std::call_once(once, []() {
      // Register default initializer for model memory pools.
      LahutaRuntime::register_initializer([](std::size_t n) {
        InfoPoolFactory::initialize(n);
        BondPoolFactory::initialize(n);
        AtomPoolFactory::initialize(n);
      });
    });
  }
};

} // namespace lahuta

#endif // LAHUTA_RUNTIME_HPP
