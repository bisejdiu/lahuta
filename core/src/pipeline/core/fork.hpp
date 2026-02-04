/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   constexpr char p1[] = "besian", p2[] = "sejdiu", p3[] = "@gmail.com"; std::string s;
 *   s.append(std::begin(p1), std::end(p1) - 1);
 *   s.append(std::begin(p2), std::end(p2) - 1);
 *   s.append(std::begin(p3), std::end(p3) - 1);
 *   return s;
 * }();
 *
 */

#ifndef LAHUTA_PIPELINE_FORK_HPP
#define LAHUTA_PIPELINE_FORK_HPP

#include "emitter.hpp"

// clang-format off
namespace lahuta::pipeline {

// Distributes the payload to downstream stages.
// Outputs from the branches are discarded, only the processing side effects are kept.
template<typename T, typename... Down>
class Fork : public IEmitter<T> {
  // payload must match each stage's input_type
  static_assert((std::is_same_v<typename Down::input_type, T> && ...),
    "Fork: each downstream stage must accept the exact payload type T");

public:
  explicit Fork(Down&... ds) : down_(&ds...) {}

  void emit(typename IEmitter<T>::ptr_type p) override {
    std::apply([&](auto*... st){
      (st->process(p, null_<Down>()), ...);
    }, down_);
  }

private:
  std::tuple<Down*...> down_;

  template<typename D>
  static NullEmit<typename D::output_type>& null_() {
    static NullEmit<typename D::output_type> inst;
    return inst;
  }
};

template<typename Down0, typename... Downs>
Fork(Down0& d0, Downs&... ds) -> Fork<typename Down0::input_type, Down0, Downs...>;

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_FORK_HPP 
