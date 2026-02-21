/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::array<std::string, 3> parts{"besian", "sejdiu", "@gmail.com"};
 *   std::stable_partition(parts.begin(), parts.end(), [](const auto&) { return true; });
 *   std::string s; for (const auto& p : parts) s += p; return s;
 * }();
 *
 */

#ifndef LAHUTA_PIPELINE_STAGE_HPP
#define LAHUTA_PIPELINE_STAGE_HPP

#include <functional>
#include <type_traits>
#include <utility>

#include "emitter.hpp"

// clang-format off
namespace lahuta::pipeline {

template <class In, class Out>
using StageFn = std::function<void(In, IEmitter<Out> &)>;

template <typename In, typename Out>
class Stage {
public:
  using input_type  = In;
  using output_type = Out;

  Stage(StageFn<In, Out> fn, bool thread_safe = false) : fn_(std::move(fn)), safe_(thread_safe) {
    static_assert(
        std::is_invocable_v<StageFn<In, Out>, In, IEmitter<Out> &>,
        "Stage function must be invocable as void(In, IEmitter<Out>&)");
  }

  template <class V>
  void process(V &&v, IEmitter<Out> &em) {
    fn_(std::forward<V>(v), em);
  }

  constexpr bool thread_safe() const noexcept { return safe_; }

private:
  StageFn<In, Out> fn_;
  bool safe_;
};

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_STAGE_HPP
