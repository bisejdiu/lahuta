#ifndef LAHUTA_CLI_TASKS_SHAPE_METRICS_SUMMARY_SINK_HPP
#define LAHUTA_CLI_TASKS_SHAPE_METRICS_SUMMARY_SINK_HPP

#include <filesystem>
#include <memory>

#include "pipeline/io/sink_iface.hpp"
#include "tasks/shape_metrics_task.hpp"

namespace lahuta::cli::shape_metrics {
namespace P = lahuta::pipeline;

[[nodiscard]] std::shared_ptr<P::IDynamicSink>
make_shape_metrics_summary_sink(std::filesystem::path output_dir,
                                std::shared_ptr<ShapeMetricsCounters> counters);

} // namespace lahuta::cli::shape_metrics

#endif // LAHUTA_CLI_TASKS_SHAPE_METRICS_SUMMARY_SINK_HPP
