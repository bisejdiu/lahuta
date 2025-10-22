#ifndef LAHUTA_CLI_COMMANDS_REPORTING_HPP
#define LAHUTA_CLI_COMMANDS_REPORTING_HPP

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
