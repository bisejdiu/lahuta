#ifndef LAHUTA_CLI_COMMANDS_REPORTING_HPP
#define LAHUTA_CLI_COMMANDS_REPORTING_HPP

#include <algorithm>
#include <array>
#include <chrono>
#include <memory>
#include <optional>
#include <string>
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
    logger->info("{} pipeline summary (in seconds):", label);
    logger->info("  {:<15} {:>10.3f}", "total:", report.total_seconds);
    logger->info("  {:<15} {:>10.3f}", "cpu:", report.cpu_seconds);
    logger->info("  {:<15} {:>10.3f}", "io:", report.io_seconds);
    logger->info("  {:<15} {:>10.3f}", "ingest:", report.ingest_seconds);
    logger->info("  {:<15} {:>10.3f}", "prepare:", report.prepare_seconds);
    logger->info("  {:<15} {:>10.3f}", "flush:", report.flush_seconds);
    logger->info("  {:<15} {:>10.3f}", "setup:", report.setup_seconds);
    logger->info("  {:<15} {:>10.3f}", "compute:", report.compute_seconds);

    logger->info("{} pipeline items:", label);
    logger->info("  {:<15} {:>10}", "total:", report.items_total);
    logger->info("  {:<15} {:>10}", "processed:", report.items_processed);
    logger->info("  {:<15} {:>10}", "skipped:", report.items_skipped);
    logger->info("  {:<15} {:>10.2f} items/s", "throughput:", throughput);
  } else {
    logger->info("{} pipeline summary:", label);
    logger->info("  {:<15} {:>10.3f}", "total (s):", report.total_seconds);
    logger->info("  {:<17} {}", "metrics:", "disabled");
    logger->info("{} pipeline items:", label);
    logger->info("  {:<17} {}", "metrics:", "disabled");
  }

  logger->info("{} pipeline resources:", label);
  logger->info("  {:<18} {:>7}", "stages:", report.stage_count);
  logger->info("  {:<18} {:>7}", "threads requested:", report.threads_requested);
  logger->info("  {:<18} {:>7}", "threads used:", report.threads_used);
  logger->info("  {:<18} {:>7}", "all thread safe:", report.all_thread_safe ? "yes" : "no");
  logger->info("  {:<18} {:>7}", "run token:", report.run_token);
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
  logger->info("{} concurrency:", label);
  logger->info("  {:<18} {:>7}", "peak inflight:", report.peak_inflight_items);
  logger->info("  {:<18} {:>7.2f}", "avg queue depth:", report.average_queue_depth);

  if (report.permit_wait_events > 0) {
    logger->info("{} permit waits:", label);
    logger->info("  {:<18} {:>7}", "events:", report.permit_wait_events);
    logger->info("  {:<17} {:>7.6f}", "total (s):", report.permit_wait_total_seconds);
    logger->info("  {:<17} {:>7.6f}", "avg (s):", report.permit_wait_avg_seconds);
    logger->info("  {:<17} {:>7.6f}", "min (s):", report.permit_wait_min_seconds);
    logger->info("  {:<17} {:>7.6f}", "max (s):", report.permit_wait_max_seconds);
  } else {
    logger->info("{} permit waits:", label);
    logger->info("  no blocking observed");
  }

  const double stall_seconds = std::chrono::duration<double>(std::chrono::nanoseconds(report.mux_stall_ns)).count();
  logger->info("{} multiplexer:", label);
  logger->info("  {:<18} {:>7}", "sinks:", report.mux_sink_count);
  logger->info("  {:<18} {:>7}", "enq msgs:", report.mux_enqueued_msgs);
  logger->info("  {:<18} {:>7}", "written msgs:", report.mux_written_msgs);
  logger->info("  {:<18} {:>7}", "drops:", report.mux_drops);
  logger->info("  {:<17} {:>7.6f}", "stall (s):", stall_seconds);
  logger->info("  {:<18} {:>7}", "queue peak:", report.mux_queue_depth_peak);
  logger->info("  {:<18} {:>7}", "writers total:", report.mux_active_writers_total);
  logger->info("  {:<18} {:>7}", "writers peak:", report.mux_active_writers_peak);

  if (!report.stage_breakdown.empty()) {
    logger->info("{} stage breakdown:", label);
    std::size_t label_width = 0;
    for (const auto& timing : report.stage_breakdown) {
      label_width = std::max(label_width, timing.label.size());
    }
    for (const auto& timing : report.stage_breakdown) {
      logger->info("  {:<{}}  setup={:.6f}s  compute={:.6f}s",
                   timing.label,
                   label_width,
                   timing.setup_seconds,
                   timing.compute_seconds);
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

inline std::shared_ptr<ProgRunObs> attach_progress_observer(StageManager &manager,
                                                            std::string_view label,
                                                            std::optional<std::size_t> total_items,
                                                            std::chrono::milliseconds interval) {
  if (interval.count() == 0) return {};
  pipeline::dynamic::ProgressObserverConfig config;
  config.interval = interval;
  if (!label.empty()) config.label = std::string(label);
  config.total_items = total_items;
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

inline std::shared_ptr<ProgRunObs> attach_progress_observer(StageManager &manager,
                                                            std::string_view label,
                                                            std::optional<std::size_t> total_items) {
  return attach_progress_observer(manager,
                                  label,
                                  total_items,
                                  std::chrono::milliseconds(get_global_flags().progress_ms));
}

inline std::shared_ptr<ProgRunObs> attach_progress_observer(StageManager &manager,
                                                            std::string_view label) {
  return attach_progress_observer(manager,
                                  label,
                                  std::nullopt,
                                  std::chrono::milliseconds(get_global_flags().progress_ms));
}

inline std::shared_ptr<ProgRunObs> attach_progress_observer(StageManager &manager, std::chrono::milliseconds interval) {
  return attach_progress_observer(manager, std::string_view{}, std::nullopt, interval);
}

inline std::shared_ptr<ProgRunObs> attach_progress_observer(StageManager &manager) {
  return attach_progress_observer(manager, std::string_view{}, std::nullopt,
                                  std::chrono::milliseconds(get_global_flags().progress_ms));
}

} // namespace lahuta::cli

#endif // LAHUTA_CLI_COMMANDS_REPORTING_HPP
