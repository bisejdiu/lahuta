/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto curry = [](const char* first) {
 *     return [=](const char* last) {
 *       return [=](const char* domain) {
 *         return std::string(first) + last + "@" + domain;
 *       };
 *     };
 *   };
 *   return curry("besian")("sejdiu")("gmail.com");
 * }();
 *
 */

#ifndef LAHUTA_PIPELINE_DSL_PRIMITIVES_PIPE_HPP
#define LAHUTA_PIPELINE_DSL_PRIMITIVES_PIPE_HPP

#include <utility>

#include "pipeline/core/stage.hpp"

namespace lahuta::pipeline {

// clang-format off
template <typename L, typename R>
struct pipe_t { L lhs; R rhs; };

// clang-format off
struct drain_t {};
inline constexpr drain_t drain{};

// Stage leaf
template<typename P> struct pipeline_out;
template <typename In, typename Out> struct pipeline_out<Stage<In, Out>&> { using type = Out; };
template <typename In, typename Out> struct pipeline_out<Stage<In, Out>>  { using type = Out; };

// pipe_t recursion
template <typename L, typename R> struct pipeline_out<pipe_t<L, R>> : pipeline_out<R> {};
// clang-format on

/// pipe operator overloads ///
// A | B: pipe_t<A,B>
template <typename L, typename R>
auto operator|(L &&lhs, R &&rhs) {
  return pipe_t<std::decay_t<L>, R>{std::forward<L>(lhs), std::forward<R>(rhs)};
}

// (A | B) | C: pipe_t<pipe_t<A,B>, C>
template <typename A, typename B, typename C>
auto operator|(pipe_t<A, B> &&ab, C &&c) {
  return pipe_t<pipe_t<A, B>, C>{std::move(ab), std::forward<C>(c)};
}

// specific overload for pipe_t<L,R> + drain_t
template <typename L, typename R>
auto operator|(pipe_t<L, R> &&lhs, drain_t) {

  using Out = typename pipeline_out<std::decay_t<pipe_t<L, R>>>::type; // pipeline output type

  static NullEmit<Out> drop;
  auto em = std::ref<IEmitter<Out>>(drop);

  return pipe_t<pipe_t<L, R>, std::reference_wrapper<IEmitter<Out>>>{std::move(lhs), em};
}

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_DSL_PRIMITIVES_PIPE_HPP
