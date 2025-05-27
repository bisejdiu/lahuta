#include "cli/global_options.hpp"
#include "cli/arg_validation.hpp"
#include "commands/run.hpp"
#include "logging.hpp"
#include <iostream>
#include <string_view>

namespace lahuta::cli {

namespace global_opts {
const option::Descriptor usage[] = {
  {GlobalOptionIndex::Unknown, 0, "", "", validate::Unknown,
   "Usage: lahuta [global-options] <subcommand> [subcommand-options]\n\n"
   "Commands:\n"
   "  run <input_file> [options]    Compute inter-atomic contacts\n\n"
   "Main Options:"},
  {GlobalOptionIndex::Help, 0, "h", "help", option::Arg::None,
   "  --help,  -h                  \tPrint this help message and exit."},
  {GlobalOptionIndex::Verbose, 0, "v", "verbose", validate::Verbosity,
   "  --verbose, -v <level>        \tSet verbosity level:\n"
   "                               \t  0 = errors only\n"
   "                               \t  1 = warnings and errors (default)\n"
   "                               \t  2 = info, warnings, errors, and debug"},
  {0, 0, 0, 0, 0, 0}
};
} // namespace global_opts

[[nodiscard]] const std::unordered_map<std::string, CommandFactory>& get_command_registry() noexcept {
  static const std::unordered_map<std::string, CommandFactory> registry = {
    {"run", &RunCommand::create}
  };
  return registry;
}

[[nodiscard]] std::string parse_global_options(int argc, char* argv[], lahuta::Logger::LogLevel& log_level, int& sub_argc, char**& sub_argv) {
  argc -= (argc > 0);
  argv += (argc > 0);

  if (argc == 0) {
    option::printUsage(std::cout, global_opts::usage);
    return {};
  }

  if (argc > 0 && (std::string_view{argv[0]} == "-h" || std::string_view{argv[0]} == "--help")) {
    option::printUsage(std::cout, global_opts::usage);
    return {};
  }

  // Default log level is Warn (level 1)
  log_level = lahuta::Logger::LogLevel::Warn;
  int subcommand_start = 0;

  // Look for global options before the subcommand
  for (int i = 0; i < argc; ++i) {
    const std::string_view arg{argv[i]};
    if (arg == "-v" || arg == "--verbose") {
      // Check if there's a next argument for the verbosity level
      if (i + 1 < argc) {
        const std::string_view level{argv[i + 1]};
        if (level == "0") {
          log_level = lahuta::Logger::LogLevel::Error;
        } else if (level == "1") {
          log_level = lahuta::Logger::LogLevel::Warn;
        } else if (level == "2") {
          log_level = lahuta::Logger::LogLevel::Debug;
        } else {
          lahuta::Logger::get_logger()->error("Invalid verbosity level '{}'. Must be 0 (errors only), 1 (warnings+), or 2 (info+debug){}", level, validate::HELP_MSG_SUFFIX);
          return {};
        }
        ++i; // Skip the verbosity level argument
      } else {
        lahuta::Logger::get_logger()->error("Option '{}' requires a verbosity level (0, 1, or 2){}", arg, validate::HELP_MSG_SUFFIX);
        return {};
      }
    } else if (arg == "-h" || arg == "--help") {
      // Help option - already handled above, but include here for completeness
      option::printUsage(std::cout, global_opts::usage);
      return {};
    } else if (arg.empty() || arg[0] != '-') {
      // Found the subcommand
      subcommand_start = i;
      break;
    } else {
      // Unknown global option
      lahuta::Logger::get_logger()->error("Unknown global option '{}'{}", arg, validate::HELP_MSG_SUFFIX);
      return {};
    }
  }

  if (subcommand_start >= argc) {
    lahuta::Logger::get_logger()->error("Subcommand is required{}", validate::HELP_MSG_SUFFIX);
    return {};
  }

  const std::string subcommand{argv[subcommand_start]}; // get subcommand

  // prep args for subcommand
  sub_argc = argc - subcommand_start - 1;
  sub_argv = argv + subcommand_start + 1;

  return subcommand;
}

} // namespace lahuta::cli
