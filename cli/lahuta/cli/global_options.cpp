#include <iostream>
#include <string_view>

#include "cli/arg_validation.hpp"
#include "cli/global_options.hpp"
#include "commands/contacts.hpp"
#include "commands/createdb.hpp"
#include "logging.hpp"

// clang-format off
namespace lahuta::cli {
namespace global_opts {

const option::Descriptor usage[] = {
  {GlobalOptionIndex::Unknown, 0, "", "", validate::Ignore,
   "Usage: lahuta [global-options] [subcommand-options]\n\n"
   "Lahuta computes inter-atomic contacts using high-performance pipeline architecture.\n"
   "The 'contacts' subcommand is used by default when no subcommand is specified.\n\n"
   "Available subcommands:\n"
   "  contacts [options]            Compute inter-atomic contacts (default)\n"
   "  createdb [options]            Create database from structure files\n\n"
   "Main Options:"},
  {GlobalOptionIndex::Help, 0, "h", "help", option::Arg::None,
   "  --help,  -h                  \tPrint this help message and exit."},
  {GlobalOptionIndex::Verbose, 0, "v", "verbose", validate::Verbosity,
   "  --verbose, -v <level>        \tSet verbosity level:\n"
   "                               \t  0 = errors only\n"
   "                               \t  1 = info, warnings, and errors\n"
   "                               \t  2 = debug, info, warnings, and errors"},
  {0, 0, 0, 0, 0, 0}
};
} // namespace global_opts

[[nodiscard]] const std::unordered_map<std::string, CommandFactory>& get_command_registry() noexcept {
  static const std::unordered_map<std::string, CommandFactory> registry = {
    {"contacts", &ContactsCommand::create},
    {"createdb", &CreateDbCommand::create}
  };
  return registry;
}

void print_global_help() {
  option::printUsage(std::cout, global_opts::usage);
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

  // Default log level is Warn unless overridden by -v
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
          log_level = lahuta::Logger::LogLevel::Info;
        } else if (level == "2") {
          log_level = lahuta::Logger::LogLevel::Debug;
        } else {
          lahuta::Logger::get_logger()->error("Invalid verbosity level '{}'. Must be 0 (errors only), 1 (info+), or 2 (debug+){}", level, validate::HELP_MSG_SUFFIX);
          return {};
        }
        ++i; // Skip the verbosity level argument
      } else {
        lahuta::Logger::get_logger()->error("Option '{}' requires a verbosity level (0, 1, or 2){}", arg, validate::HELP_MSG_SUFFIX);
        return {};
      }
    } else if (arg == "-h" || arg == "--help") {
      // Help option, already handled above
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
    // No subcommand provided then default to "contacts"
    sub_argc = argc;
    sub_argv = argv;
    return "contacts";
  }

  const std::string subcommand{argv[subcommand_start]}; // get subcommand

  // Check if this is actually a subcommand or just an option for the default contacts command
  const auto& registry = get_command_registry();
  if (registry.find(subcommand) == registry.end()) {
    // TODO: re-evaluate if we want this behavior.
    // Not a recognized subcommand, treat as options for default "contacts" command
    sub_argc = argc;
    sub_argv = argv;
    return "contacts";
  }

  // prep args for subcommand
  sub_argc = argc - subcommand_start - 1;
  sub_argv = argv + subcommand_start + 1;

  return subcommand;
}

} // namespace lahuta::cli
