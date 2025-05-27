#include "cli/global_options.hpp"
#include "cli/arg_validation.hpp"
#include "commands/run.hpp"
#include "logging.hpp"
#include <iostream>

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
  {GlobalOptionIndex::Quiet, 0, "q", "quiet", option::Arg::None,
   "  --quiet, -q                  \tSuppress verbose output."},
  {0, 0, 0, 0, 0, 0}
};
} // namespace global_opts

[[nodiscard]] const std::unordered_map<std::string, CommandFactory>& get_command_registry() noexcept {
  static const std::unordered_map<std::string, CommandFactory> registry = {
    {"run", &RunCommand::create}
  };
  return registry;
}

[[nodiscard]] std::string parse_global_options(int argc, char* argv[], bool& global_quiet, int& sub_argc, char**& sub_argv) {
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

  global_quiet = false;
  int subcommand_start = 0;

  // Look for global options before the subcommand
  for (int i = 0; i < argc; ++i) {
    const std::string_view arg{argv[i]};
    if (arg == "-q" || arg == "--quiet") {
      global_quiet = true;
    } else if (arg.empty() || arg[0] != '-') {
      subcommand_start = i;
      break;
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
