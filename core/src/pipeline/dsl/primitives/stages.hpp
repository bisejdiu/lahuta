#ifndef LAHUTA_PIPELINE_DSL_PRIMITIVES_STAGES_HPP
#define LAHUTA_PIPELINE_DSL_PRIMITIVES_STAGES_HPP

#include <type_traits>

#include "pipeline/core/stage.hpp"
#include "pipeline/dsl/function_traits.hpp"
#include "tags.hpp"

namespace lahuta::pipeline {

// clang-format off
template<typename U> struct emitter_payload;
template<typename T> struct emitter_payload<IEmitter<T>&> { using type = T; };

template<typename Fn> using stage_in_t  = std::remove_reference_t <util::function_arg_t<Fn, 0>>;
template<typename Fn> using stage_out_t = typename emitter_payload<util::function_arg_t<Fn, 1>>::type;
// clang-format on

/// stage(tag, lambda), auto deduce In/Out
namespace detail {
template <typename Tag, typename Fn>
auto make_stage(Tag, Fn &&fn) {
  using In  = stage_in_t<Fn>;
  using Out = stage_out_t<Fn>;

  static_assert(std::is_invocable_v<Fn, In, IEmitter<Out> &>, "lambda must be void(In, IEmitter<Out>&)");

  return Stage<In, Out>{std::forward<Fn>(fn), Tag::value};
}
} // namespace detail

template <typename Fn>
auto stage(Fn &&fn) {
  return detail::make_stage(thread_unsafe, std::forward<Fn>(fn));
}
template <typename Fn>
auto stage(thread_safe_t, Fn &&fn) {
  return detail::make_stage(thread_safe, std::forward<Fn>(fn));
}
template <typename Fn>
auto stage(thread_unsafe_t, Fn &&fn) {
  return detail::make_stage(thread_unsafe, std::forward<Fn>(fn));
}

template <typename In, typename Out, typename Fn, typename Tag = thread_unsafe_t>
auto stage_typed(Fn &&fn, Tag tag = {}) {
  static_assert(std::is_invocable_v<Fn, In, IEmitter<Out> &>);
  return Stage<In, Out>{std::forward<Fn>(fn), Tag::value};
}

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_DSL_PRIMITIVES_STAGES_HPP
