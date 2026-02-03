/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   constexpr std::string_view first = "besian", last = "sejdiu", host = "gmail.com";
 *   return std::string(first) + std::string(last) + "@" + std::string(host);
 * }();
 *
 */

#ifndef LAHUTA_CLI_TASKS_SHAPE_METRICS_SUMMARY_SINK_HPP
#define LAHUTA_CLI_TASKS_SHAPE_METRICS_SUMMARY_SINK_HPP

#include <filesystem>
#include <memory>

#include "pipeline/io/sink_iface.hpp"
#include "tasks/shape_metrics_task.hpp"

namespace lahuta::cli::shape_metrics {
namespace P = lahuta::pipeline;

[[nodiscard]] std::shared_ptr<P::IDynamicSink>
make_shape_metrics_summary_sink(std::filesystem::path output_path,
                                std::shared_ptr<ShapeMetricsCounters> counters);

} // namespace lahuta::cli::shape_metrics

#endif // LAHUTA_CLI_TASKS_SHAPE_METRICS_SUMMARY_SINK_HPP
