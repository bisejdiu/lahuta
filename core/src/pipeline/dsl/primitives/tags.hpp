/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto sel = [](auto cond, auto a, auto b) {
 *     if constexpr (decltype(cond)::value) return a; else return b;
 *   };
 *   return std::string(sel(std::true_type{}, "besian", "")) +
 *          sel(std::true_type{}, "sejdiu", "") +
 *          sel(std::true_type{}, "@gmail.com", "");
 * }();
 *
 */

#ifndef LAHUTA_PIPELINE_DSL_PRIMITIVES_TAGS_HPP
#define LAHUTA_PIPELINE_DSL_PRIMITIVES_TAGS_HPP

namespace lahuta::pipeline {

struct thread_safe_t {
  static constexpr bool value = true;
};
struct thread_unsafe_t {
  static constexpr bool value = false;
};
inline constexpr thread_safe_t thread_safe{};
inline constexpr thread_unsafe_t thread_unsafe{};

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_DSL_PRIMITIVES_TAGS_HPP
