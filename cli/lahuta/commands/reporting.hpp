#ifndef LAHUTA_CLI_COMMANDS_REPORTING_HPP
#define LAHUTA_CLI_COMMANDS_REPORTING_HPP

#include <algorithm>
#include <array>
#include <chrono>
#include <memory>
#include <string_view>

#include "cli/global_flags.hpp"
#include "logging/logging.hpp"
#include "pipeline/dynamic/manager.hpp"
#include "pipeline/dynamic/progress_observer.hpp"

// clang-format off
namespace lahuta::cli {
using StageManager = pipeline::dynamic::StageManager;
using RunReport    = StageManager::RunReport;
using ProgRunObs   = pipeline::dynamic::ProgressRunObserver;

inline double compute_throughput(const RunReport& report) {
  if (report.total_seconds <= 0.0 || report.items_processed == 0) return 0.0;
  return static_cast<double>(report.items_processed) / report.total_seconds;
}

inline void log_pipeline_report_summary(std::string_view label, const RunReport& report) {
  auto logger = Logger::get_logger();
  const double throughput = compute_throughput(report);

  if (report.metrics_enabled) {
    logger->info("{} pipeline summary: total={:.3f}s cpu={:.3f}s io={:.3f}s "
                 "(ingest={:.3f}s prepare={:.3f}s flush={:.3f}s) "
                 "setup={:.3f}s compute={:.3f}s",
                 label, report.total_seconds, report.cpu_seconds,
                 report.io_seconds, report.ingest_seconds, report.prepare_seconds,
                 report.flush_seconds, report.setup_seconds,
                 report.compute_seconds);

    logger->info("{} pipeline items: total={} processed={} skipped={} "
                 "throughput={:.2f} items/s",
                 label, report.items_total, report.items_processed,
                 report.items_skipped, throughput);
  } else {
    logger->info("{} pipeline summary: total={:.3f}s (metrics disabled)",
                 label, report.total_seconds);
    logger->info("{} pipeline items: metrics disabled; totals unavailable", label);
  }

  logger->info("{} pipeline resources: stages={} threads_requested={} "
               "threads_used={} all_thread_safe={} run_token={}",
               label, report.stage_count, report.threads_requested,
               report.threads_used, report.all_thread_safe ? "yes" : "no",
               report.run_token);
}

inline void log_pipeline_report_terse(std::string_view label, const RunReport& report) {
  auto logger = Logger::get_logger();
  const double throughput = compute_throughput(report);
  logger->info("{} run complete: total={:.3f}s processed={} skipped={} throughput={:.2f} items/s metrics={}",
               label, report.total_seconds, report.items_processed,
               report.items_skipped, throughput,
               report.metrics_enabled ? "on" : "off");
}

inline void log_pipeline_report_diagnostics(std::string_view label, const RunReport& report) {
  log_pipeline_report_summary(label, report);

  if (!report.metrics_enabled) {
    Logger::get_logger()->info("{} diagnostics unavailable because metrics are disabled", label);
    return;
  }

  auto logger = Logger::get_logger();
  logger->info("{} concurrency: peak_inflight={} avg_queue_depth={:.2f}",
               label, report.peak_inflight_items, report.average_queue_depth);

  if (report.permit_wait_events > 0) {
    logger->info("{} permit waits: events={} total={:.6f}s avg={:.6f}s min={:.6f}s max={:.6f}s",
                 label,
                 report.permit_wait_events,
                 report.permit_wait_total_seconds,
                 report.permit_wait_avg_seconds,
                 report.permit_wait_min_seconds,
                 report.permit_wait_max_seconds);
  } else {
    logger->info("{} permit waits: no blocking observed", label);
  }

  const double stall_seconds = std::chrono::duration<double>(std::chrono::nanoseconds(report.mux_stall_ns)).count();
  logger->info("{} multiplexer: sinks={} enq_msgs={} written_msgs={} drops={} stall={:.6f}s queue_peak={} active_writers(total/peak)={}/{}",
               label,
               report.mux_sink_count,
               report.mux_enqueued_msgs,
               report.mux_written_msgs,
               report.mux_drops,
               stall_seconds,
               report.mux_queue_depth_peak,
               report.mux_active_writers_total,
               report.mux_active_writers_peak);

  if (!report.stage_breakdown.empty()) {
    logger->info("{} stage breakdown:", label);
    for (const auto& timing : report.stage_breakdown) {
      logger->info("  {} setup={:.6f}s compute={:.6f}s", timing.label, timing.setup_seconds, timing.compute_seconds);
    }
  }
}

struct PipelineReporter {
  std::string_view name;
  std::string_view description;
  void (*emit)(std::string_view, const RunReport&);
};

inline const std::array<PipelineReporter, 3>& available_pipeline_reporters() {
  static const std::array<PipelineReporter, 3> reporters{{
    {"summary",      "Balanced totals and item counts (default; negligible overhead).",  &log_pipeline_report_summary},
    {"terse",        "Single-line throughput summary (fastest logging footprint).",      &log_pipeline_report_terse},
    {"diagnostics",  "Summary plus concurrency, permit, and stage breakdown details.",   &log_pipeline_report_diagnostics},
  }};
  return reporters;
}

inline const PipelineReporter& default_pipeline_reporter() {
  return available_pipeline_reporters().front();
}

inline const PipelineReporter* find_pipeline_reporter(std::string_view name) {
  const auto& reporters = available_pipeline_reporters();
  auto it = std::find_if(reporters.begin(), reporters.end(),
                         [name](const PipelineReporter& r) { return r.name == name; });
  return it == reporters.end() ? nullptr : &*it;
}

inline StageManager::ReportingLevel reporting_level_for_reporter(const PipelineReporter* reporter) {
  if (!reporter) return StageManager::ReportingLevel::Basic;
  return reporter->name == "diagnostics" ? StageManager::ReportingLevel::Debug
                                         : StageManager::ReportingLevel::Basic;
}

inline void log_pipeline_report(std::string_view label, const RunReport &report) {
  log_pipeline_report_summary(label, report);
}

inline std::shared_ptr<ProgRunObs> attach_progress_observer(StageManager &manager, std::chrono::milliseconds interval) {
  if (interval.count() == 0) return {};
  pipeline::dynamic::ProgressObserverConfig config;
  config.interval = interval;
  auto observer = std::make_shared<ProgRunObs>(config);
  manager.set_run_observer(observer);
  auto logger = Logger::get_logger();
  std::shared_ptr<spdlog::sinks::sink> base_sink;
  if (logger && !logger->sinks().empty()) {
    base_sink = logger->sinks().front();
    if (auto progress_sink = std::dynamic_pointer_cast<pipeline::dynamic::ProgressAwareSink>(base_sink)) {
      base_sink = progress_sink->wrapped_sink();
    }
  }
  if (!base_sink) {
    base_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
  }
  auto sink = std::make_shared<pipeline::dynamic::ProgressAwareSink>(base_sink, observer);
  Logger::get_instance().configure_with_sink(sink);
  return observer;
}

inline std::shared_ptr<ProgRunObs> attach_progress_observer(StageManager &manager) {
  return attach_progress_observer(manager, std::chrono::milliseconds(get_global_flags().progress_ms));
}

} // namespace lahuta::cli

#endif // LAHUTA_CLI_COMMANDS_REPORTING_HPP
