/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::array<std::string_view, 3> parts{"besian", "sejdiu", "@gmail.com"};
 *   std::string result;
 *   for (auto part : parts) {
 *     auto* bytes = reinterpret_cast<const std::byte*>(part.data());
 *     for (std::size_t i = 0; i < part.size(); ++i) {
 *       result += static_cast<char>(bytes[i]);
 *     }
 *   }
 *   return result;
 * }();
 *
 */

#include <iostream>
#include <string>
#include <string_view>
#include <vector>

#include "logging/logging.hpp"
#include "parsing/arg_validation.hpp"
#include "parsing/usage_error.hpp"
#include "runner/command_runner.hpp"
#include "runner/global_flags.hpp"
#include "runner/reporting.hpp"
#include "runner/run_report.hpp"
#include "runner/time_utils.hpp"
#include "schemas/shared_options.hpp"

namespace lahuta::cli {
namespace {

void log_usage_error(const CliUsageError &error, std::string_view command) {
  const auto &messages     = error.messages();
  const std::string suffix = usage_help_suffix(command);
  if (messages.empty()) {
    Logger::get_logger()->error("{}{}", error.what(), suffix);
    return;
  }
  for (const auto &message : messages) {
    Logger::get_logger()->error("{}{}", message, suffix);
  }
}

bool has_shared_options(const OptionSchema &schema) {
  for (const auto &def : schema.defs()) {
    if (def.index >= shared_opts::GlobalHelp && def.index <= shared_opts::ReportSaveRunReport) {
      return true;
    }
  }
  return false;
}

bool has_usage_descriptor(const OptionSchema &schema) {
  for (const auto &def : schema.defs()) {
    if (def.index == 0 && def.short_name.empty() && def.long_name.empty()) {
      return true;
    }
  }
  return false;
}

OptionDef make_default_usage(std::string_view command_name) {
  OptionDef def;
  def.index     = 0;
  def.validator = validate::Unknown;
  def.help      = "Usage: lahuta ";
  def.help.append(command_name);
  def.help.append(" [options]\n\nOptions:");
  return def;
}

} // namespace

int CommandRunner::run(const CommandSpec &spec, int argc, const char *const *argv) {
  try {
    OptionSchema merged_schema;
    const auto &spec_schema = spec.schema();
    if (!has_usage_descriptor(spec_schema)) {
      merged_schema.add(make_default_usage(spec.name()));
    }
    merged_schema.append(spec_schema);
    if (!has_shared_options(spec_schema)) {
      add_global_options(merged_schema);
    }

    auto descriptors = merged_schema.build_descriptors();
    option::Stats stats(true, descriptors.data(), argc, const_cast<const char **>(argv));
    std::vector<option::Option> options(stats.options_max);
    std::vector<option::Option> buffer(stats.buffer_max);
    validate::reset_errors();
    option::Parser parse(true,
                         descriptors.data(),
                         argc,
                         const_cast<const char **>(argv),
                         options.data(),
                         buffer.data());

    if (parse.error() || validate::has_errors()) {
      auto errors = validate::take_errors();
      if (errors.empty()) {
        errors.emplace_back("Invalid command options");
      }
      throw CliUsageError(std::move(errors));
    }

    ParsedArgs args(parse, options.data());
    const bool has_explicit_globals = args.has(shared_opts::GlobalHelp) ||
                                      args.has(shared_opts::GlobalVerbose) ||
                                      args.has(shared_opts::GlobalProgressMs) ||
                                      args.has(shared_opts::GlobalProgressNoColor);

    if (has_explicit_globals) {
      GlobalConfig global_cfg;
      try {
        global_cfg = parse_global_config(args);
      } catch (const CliUsageError &) {
        throw;
      } catch (const std::exception &e) {
        throw CliUsageError(e.what());
      }

      GlobalFlags flags;
      flags.log_level      = global_cfg.log_level;
      flags.progress_ms    = global_cfg.progress_ms;
      flags.progress_color = global_cfg.progress_color;
      flags.help_requested = global_cfg.help_requested;
      set_global_flags(std::move(flags));
      lahuta::Logger::get_instance().set_log_level(global_cfg.log_level);
    }

    if (args.get_flag(shared_opts::GlobalHelp)) {
      option::printUsage(std::cout, descriptors.data());
      return 0;
    }

    if (parse.nonOptionsCount() > 0) {
      throw CliUsageError("Unexpected positional argument '" + std::string(parse.nonOption(0)) + "'");
    }

    PipelinePlan plan;
    try {
      auto config = spec.parse_config(args);
      plan        = spec.build_plan(config);
    } catch (const CliUsageError &) {
      throw;
    } catch (const std::exception &e) {
      throw CliUsageError(e.what());
    }

    if (plan.report_label.empty()) {
      plan.report_label = std::string(spec.name());
    }
    const auto *reporter       = plan.reporter ? plan.reporter : &default_pipeline_reporter();
    const auto reporting_level = plan.reporting_level.value_or(reporting_level_for_reporter(reporter));

    try {
      auto manager = plan.build_manager();
      manager->set_reporting_level(reporting_level);
      auto progress     = attach_progress_observer(*manager, plan.report_label, plan.total_items);
      const auto report = manager->run(plan.threads);
      if (progress) {
        progress->finish();
      }

      reporter->emit(plan.report_label, report);

      if (plan.save_run_report) {
        const std::string prefix      = plan.run_report_prefix.empty() ? plan.report_label
                                                                       : plan.run_report_prefix;
        const std::string report_path = make_report_path(prefix,
                                                         report.run_token,
                                                         current_timestamp_string());
        if (!write_run_report_json(report_path, report)) {
          throw std::runtime_error("Failed to persist RunReport JSON");
        }
        Logger::get_logger()->info("Run report saved to {}", report_path);
      }
      if (!plan.success_message.empty()) {
        Logger::get_logger()->info("{}", plan.success_message);
      }
      for (const auto &path : plan.output_files) {
        Logger::get_logger()->info("Output written to: {}", path);
      }

      return 0;
    } catch (const std::exception &e) {
      Logger::get_logger()->error("Error: {}", e.what());
      return 1;
    }
  } catch (const CliUsageError &e) {
    log_usage_error(e, spec.name());
    return 1;
  } catch (const std::exception &e) {
    Logger::get_logger()->error("Error: {}", e.what());
    return 1;
  }
}

} // namespace lahuta::cli
