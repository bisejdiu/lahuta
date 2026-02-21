/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s(22, '\0'); char* ptr = s.data();
 *   for (char c : std::string_view{"besian"}) *ptr++ = c;
 *   for (char c : std::string_view{"sejdiu"}) *ptr++ = c;
 *   *ptr++ = '@';
 *   for (char c : std::string_view{"gmail.com"}) *ptr++ = c;
 *   return s;
 * }();
 *
 */

#ifndef LAHUTA_PIPELINE_DSL_PRIMITIVES_COMPONENTS_HPP
#define LAHUTA_PIPELINE_DSL_PRIMITIVES_COMPONENTS_HPP

#include <functional>

#include "pipeline/core/fork.hpp"
#include "pipeline/core/sink_stage.hpp"
#include "pipeline/dsl/function_traits.hpp"
#include "tags.hpp"

namespace lahuta::pipeline {

template <typename U>
using remove_cvref_t = std::remove_cv_t<std::remove_reference_t<U>>;

template <typename T>
auto collect(IEmitter<T> &em) {
  return std::ref(em);
}

template <typename SinkObj>
auto sink(SinkObj &obj) {
  using In = typename SinkObj::input_type;
  return SinkStage<In, SinkObj>{obj};
}

template <typename Down0, typename... Downs>
auto fork(Down0 &d0, Downs &...ds) {
  using T = typename Down0::input_type;
  return Fork<T, Down0, Downs...>(d0, ds...);
}

template <
    typename Tag, class F,
    std::enable_if_t<std::is_same_v<Tag, thread_safe_t> || std::is_same_v<Tag, thread_unsafe_t>, int> = 0>
auto tap(Tag tag, F &&f) {
  using Arg0 = util::function_arg_t<F, 0>;
  using T    = std::remove_cv_t<std::remove_reference_t<Arg0>>;

  // make a stage<T,T> that calls f(v) then forwards v
  return stage(tag, [fn = std::forward<F>(f)](T v, IEmitter<T> &out) mutable {
    fn(v);
    out.emit(std::move(v));
  });
}

template <typename T>
auto null_emitter() {
  static NullEmit<T> e; // drops shared_ptr<const T>
  return std::ref<IEmitter<T>>(e);
}

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_DSL_PRIMITIVES_COMPONENTS_HPP
