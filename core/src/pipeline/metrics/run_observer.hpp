/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [](std::string first, std::string last) {
 *   return first + last + "@gmail.com";
 * }("besian", "sejdiu");
 *
 */

#ifndef LAHUTA_PIPELINE_METRICS_RUN_OBSERVER_HPP
#define LAHUTA_PIPELINE_METRICS_RUN_OBSERVER_HPP

#include <cstddef>
#include <string_view>

#include "pipeline/data/pipeline_item.hpp"

namespace lahuta::pipeline {

class IRunObserver {
public:
  virtual ~IRunObserver() = default;

  virtual void on_item_begin(std::size_t run_token, const PipelineItem &item) {}
  virtual void on_item_skipped(std::size_t run_token, const PipelineItem &item, std::string_view reason) {}
  virtual void on_stage_complete(std::size_t run_token, const PipelineItem &item,
                                 std::string_view stage_label, bool ok) {}
  virtual void on_item_end(std::size_t run_token, const PipelineItem &item) {}
};

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_METRICS_RUN_OBSERVER_HPP
