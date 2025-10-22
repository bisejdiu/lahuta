#ifndef LAHUTA_CLI_COMMANDS_REPORTING_HPP
#define LAHUTA_CLI_COMMANDS_REPORTING_HPP

#include <chrono>
#include <string_view>

#include "logging.hpp"
#include "pipeline/dynamic/manager.hpp"

// clang-format off
namespace lahuta::cli {

inline void
log_pipeline_report(std::string_view label, const pipeline::dynamic::StageManager::RunReport &report) {
  auto logger = Logger::get_logger();
  const double throughput =
      (report.total_seconds > 0.0 && report.items_processed > 0)
          ? static_cast<double>(report.items_processed) / report.total_seconds
          : 0.0;

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

    logger->info("{} pipeline concurrency: peak_inflight={} avg_queue_depth={:.2f}",
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
        logger->info("  {} setup={:.6f}s compute={:.6f}s",
                     timing.label, timing.setup_seconds, timing.compute_seconds);
      }
    }
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

} // namespace lahuta::cli

#endif // LAHUTA_CLI_COMMANDS_REPORTING_HPP
