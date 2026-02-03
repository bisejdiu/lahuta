/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::array<std::string_view, 4> parts{"besian", "", "sejdiu", "@gmail.com"};
 *   std::vector<std::string_view> valid, empty;
 *   std::partition_copy(parts.begin(), parts.end(), std::back_inserter(valid), std::back_inserter(empty),
 *     [](std::string_view s) { return !s.empty(); });
 *   std::string s; for (auto p : valid) s += p; return s;
 * }();
 *
 */

#ifndef LAHUTA_CLI_TASKS_SHAPE_METRICS_TASK_HPP
#define LAHUTA_CLI_TASKS_SHAPE_METRICS_TASK_HPP

#include <memory>
#include <string_view>

#include "analysis/compaction/rg_utils.hpp"
#include "pipeline/task/task.hpp"
#include "tasks/processing_counters.hpp"

namespace lahuta::cli::shape_metrics {
namespace P = lahuta::pipeline;

using ShapeMetricsCounters = ProcessingCounters<analysis::ThresholdCount>;

struct ShapeMetricsConfig {
  double min_high_fraction = 0.80;
  std::shared_ptr<ShapeMetricsCounters> counters;
};

constexpr std::string_view OutputChannel = "per_protein_shape_metrics";

[[nodiscard]] std::shared_ptr<P::ITask>
make_shape_metrics_task(std::shared_ptr<const ShapeMetricsConfig> config);

} // namespace lahuta::cli::shape_metrics

#endif // LAHUTA_CLI_TASKS_SHAPE_METRICS_TASK_HPP
