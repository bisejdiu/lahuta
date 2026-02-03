/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   constexpr std::array parts{"besian", "sejdiu", "@gmail.com"};
 *   return std::accumulate(parts.begin(), parts.end(), std::string{});
 * }();
 *
 */

#ifndef LAHUTA_CLI_RUN_REPORT_HPP
#define LAHUTA_CLI_RUN_REPORT_HPP

#include <fstream>
#include <string>

#include "logging/logging.hpp"
#include "pipeline/runtime/manager.hpp"
#include "serialization/json.hpp"

namespace lahuta::cli {
namespace P = lahuta::pipeline;

inline std::string make_report_path(const std::string &prefix, std::size_t run_token,
                                    const std::string &timestamp) {
  return prefix + "_run_report_" + std::to_string(run_token) + '_' + timestamp + ".json";
}

inline bool write_run_report_json(const std::string &path, const P::StageManager::RunReport &report) {
  std::ofstream out(path, std::ios::out | std::ios::trunc);
  if (!out) {
    Logger::get_logger()->error("Unable to write RunReport JSON to '{}'", path);
    return false;
  }

  JsonBuilder json(512);
  json.key("total_seconds")
      .value(report.total_seconds)
      .key("cpu_seconds")
      .value(report.cpu_seconds)
      .key("io_seconds")
      .value(report.io_seconds)
      .key("ingest_seconds")
      .value(report.ingest_seconds)
      .key("prepare_seconds")
      .value(report.prepare_seconds)
      .key("flush_seconds")
      .value(report.flush_seconds)
      .key("setup_seconds")
      .value(report.setup_seconds)
      .key("compute_seconds")
      .value(report.compute_seconds)
      .key("items_total")
      .value(report.items_total)
      .key("items_processed")
      .value(report.items_processed)
      .key("items_skipped")
      .value(report.items_skipped)
      .key("stage_count")
      .value(report.stage_count)
      .key("threads_requested")
      .value(report.threads_requested)
      .key("threads_used")
      .value(report.threads_used)
      .key("all_thread_safe")
      .value(report.all_thread_safe)
      .key("run_token")
      .value(report.run_token)
      .key("metrics_enabled")
      .value(report.metrics_enabled)
      .key("peak_inflight_items")
      .value(report.peak_inflight_items)
      .key("average_queue_depth")
      .value(report.average_queue_depth)
      .key("permit_wait_total_seconds")
      .value(report.permit_wait_total_seconds)
      .key("permit_wait_min_seconds")
      .value(report.permit_wait_min_seconds)
      .key("permit_wait_max_seconds")
      .value(report.permit_wait_max_seconds)
      .key("permit_wait_avg_seconds")
      .value(report.permit_wait_avg_seconds)
      .key("permit_wait_events")
      .value(report.permit_wait_events)
      .key("mux_sink_count")
      .value(report.mux_sink_count)
      .key("mux_enqueued_msgs")
      .value(report.mux_enqueued_msgs)
      .key("mux_enqueued_bytes")
      .value(report.mux_enqueued_bytes)
      .key("mux_written_msgs")
      .value(report.mux_written_msgs)
      .key("mux_written_bytes")
      .value(report.mux_written_bytes)
      .key("mux_stall_ns")
      .value(report.mux_stall_ns)
      .key("mux_drops")
      .value(report.mux_drops)
      .key("mux_queue_depth_peak")
      .value(report.mux_queue_depth_peak)
      .key("mux_queue_bytes_peak")
      .value(report.mux_queue_bytes_peak)
      .key("mux_active_writers_total")
      .value(report.mux_active_writers_total)
      .key("mux_active_writers_peak")
      .value(report.mux_active_writers_peak);

  json.key("stage_breakdown").begin_array();
  for (const auto &stage : report.stage_breakdown) {
    json.begin_object()
        .key("label")
        .value(stage.label)
        .key("setup_seconds")
        .value(stage.setup_seconds)
        .key("compute_seconds")
        .value(stage.compute_seconds)
        .end_object();
  }
  json.end_array();

  out << json.str() << '\n';
  return true;
}

} // namespace lahuta::cli

#endif // LAHUTA_CLI_RUN_REPORT_HPP
