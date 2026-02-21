/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::array<std::string_view, 3> parts{"@gmail.com", "besian", "sejdiu"};
 *   std::sort(parts.begin(), parts.end());
 *   return std::string(parts[1]) + std::string(parts[2]) + std::string(parts[0]);
 * }();
 *
 */

#include <iostream>
#include <string>
#include <string_view>
#include <vector>

#include "logging/logging.hpp"
#include "parsing/arg_validation.hpp"
#include "parsing/global_args.hpp"
#include "parsing/parsed_args.hpp"
#include "runner/app.hpp"
#include "runner/command_runner.hpp"
#include "runner/global_flags.hpp"
#include "schemas/option_schema.hpp"
#include "schemas/shared_options.hpp"
#include "specs/command_spec.hpp"
#include "version.hpp"

namespace {

void log_usage_error(const lahuta::cli::CliUsageError &error, std::string_view command = {}) {
  const auto &messages     = error.messages();
  const std::string suffix = lahuta::cli::usage_help_suffix(command);
  if (messages.empty()) {
    lahuta::Logger::get_logger()->error("{}{}", error.what(), suffix);
    return;
  }
  for (const auto &message : messages) {
    lahuta::Logger::get_logger()->error("{}{}", message, suffix);
  }
}

} // namespace

namespace lahuta::cli {

int run(int argc, char *argv[]) {
  try {
    const auto &registry = get_command_specs();
    const auto split     = split_global_args(argc, argv);

    OptionSchema global_schema;
    global_schema.add({0, "", "", validate::Unknown, build_global_usage(registry)});
    add_global_options(global_schema);

    auto descriptors      = global_schema.build_descriptors();
    const int global_argc = static_cast<int>(split.global_args.size());
    option::Stats stats(true,
                        descriptors.data(),
                        global_argc,
                        const_cast<const char **>(split.global_args.data()));
    std::vector<option::Option> options(stats.options_max);
    std::vector<option::Option> buffer(stats.buffer_max);
    validate::reset_errors();
    option::Parser parse(true,
                         descriptors.data(),
                         global_argc,
                         const_cast<const char **>(split.global_args.data()),
                         options.data(),
                         buffer.data());

    if (parse.error() || validate::has_errors()) {
      auto errors = validate::take_errors();
      if (errors.empty()) {
        errors.emplace_back("Invalid global options");
      }
      throw CliUsageError(std::move(errors));
    }

    ParsedArgs args(parse, options.data());
    const GlobalConfig global_config = parse_global_config(args);

    GlobalFlags flags;
    flags.log_level      = global_config.log_level;
    flags.progress_ms    = global_config.progress_ms;
    flags.progress_color = global_config.progress_color;
    flags.help_requested = global_config.help_requested;
    set_global_flags(std::move(flags));
    Logger::get_instance().set_log_level(global_config.log_level);

    if (global_config.help_requested) {
      option::printUsage(std::cout, descriptors.data());
      return 0;
    }

    if (global_config.version_requested) {
      std::cout << "lahuta " << lahuta::version << '\n';
      return 0;
    }

    if (split.tail.empty()) {
      Logger::get_logger()->error("No subcommand provided{}", usage_help_suffix());
      return 1;
    }

    const std::string_view subcommand{split.tail[0]};

    if (subcommand.empty() || subcommand.front() == '-') {
      Logger::get_logger()->error("Expected subcommand before '{}'{}", subcommand, usage_help_suffix());
      return 1;
    }

    const auto it = registry.find(std::string{subcommand});

    if (it == registry.end()) {
      Logger::get_logger()->error("Unknown subcommand '{}'{}", subcommand, usage_help_suffix());
      return 1;
    }

    // subcommand args
    const int sub_argc          = static_cast<int>(split.tail.size() - 1);
    const char *const *sub_argv = split.tail.data() + 1;
    return CommandRunner::run(*it->second, sub_argc, sub_argv);

  } catch (const CliUsageError &e) {
    log_usage_error(e);
    return 1;
  } catch (const std::exception &e) {
    Logger::get_logger()->error("{}{}", e.what(), usage_help_suffix());
    return 1;
  }
}

} // namespace lahuta::cli
