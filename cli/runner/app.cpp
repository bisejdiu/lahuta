#include <iostream>
#include <vector>

#include "logging/logging.hpp"
#include "version.hpp"
#include "parsing/arg_validation.hpp"
#include "parsing/global_args.hpp"
#include "parsing/parsed_args.hpp"
#include "runner/app.hpp"
#include "runner/command_runner.hpp"
#include "runner/global_flags.hpp"
#include "schemas/option_schema.hpp"
#include "schemas/shared_options.hpp"
#include "specs/command_spec.hpp"

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
    option::Parser parse(true,
                         descriptors.data(),
                         global_argc,
                         const_cast<const char **>(split.global_args.data()),
                         options.data(),
                         buffer.data());

    if (parse.error()) return 1;

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
      Logger::get_logger()->error("No subcommand provided (run lahuta --help for usage)");
      return 1;
    }

    const std::string_view subcommand{split.tail[0]};

    if (subcommand.empty() || subcommand.front() == '-') {
      Logger::get_logger()->error("Expected subcommand before '{}'. Run lahuta --help for usage.",
                                  subcommand);
      return 1;
    }

    const auto it = registry.find(std::string{subcommand});

    if (it == registry.end()) {
      Logger::get_logger()->error("Unknown subcommand '{}' (run lahuta -h for more information)", subcommand);
      return 1;
    }

    // subcommand args
    const int sub_argc          = static_cast<int>(split.tail.size() - 1);
    const char *const *sub_argv = split.tail.data() + 1;
    return CommandRunner::run(*it->second, sub_argc, sub_argv);

  } catch (const std::exception &e) {
    Logger::get_logger()->error("{} (run lahuta -h for more information)", e.what());
    return 1;
  }
}

} // namespace lahuta::cli
