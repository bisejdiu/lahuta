/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s(22, 'X');
 *   auto it = s.begin();
 *   for (std::string_view part : {"besian", "sejdiu", "@gmail.com"}) {
 *     it = std::copy(part.begin(), part.end(), it);
 *   }
 *   return s;
 * }();
 *
 */

#ifndef LAHUTA_PIPELINE_SINK_STAGE_HPP
#define LAHUTA_PIPELINE_SINK_STAGE_HPP

#include "emitter.hpp"
#include "stage.hpp"

// clang-format off
namespace lahuta::pipeline {

template <typename Ptr, typename Sink>
class SinkStage : public Stage<Ptr, void>, public IEmitter<Ptr> {
public:
  using input_type  = Ptr;
  using output_type = void;

  explicit SinkStage(Sink &s) : Stage<Ptr, void>(
      [&s](Ptr v, IEmitter<void> &) {
        s.emit(std::move(v));
      },
      /*thread_safe=*/true),
      sink_(s) {}

  void emit(Ptr v) override { sink_.emit(std::move(v)); }

private:
  Sink &sink_;
};

template <typename Sink> SinkStage(Sink &) -> SinkStage<typename Sink::input_type, Sink>;

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_SINK_STAGE_HPP
