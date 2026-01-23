#include <iostream>

#include "cli/arg_validation.hpp"
#include "cli/global_options.hpp"
#include "commands/contacts.hpp"
#include "commands/createdb.hpp"
#include "commands/extract.hpp"
#include "commands/positions.hpp"

// clang-format off
namespace lahuta::cli {
namespace global_opts {

const option::Descriptor usage[] = {
  {GlobalOptionIndex::Unknown, 0, "", "", validate::Ignore,
   "Usage: lahuta [global-options] <subcommand> [subcommand-options]\n\n"
   "Lahuta computes inter-atomic contacts using different build-in contact providers.\n\n"
   "Available subcommands:\n"
   "  contacts  [options]            Compute inter-atomic contacts\n"
   "  createdb  [options]            Create database from structure files\n"
   "  extract   [options]            Extract data from AlphaFold2 models\n"
   "  positions [options]            Extract 3D atomic coordinates to sharded binary files\n\n"
   "Main Options:"},
  {GlobalOptionIndex::Help, 0, "h", "help", option::Arg::None,
   "  --help,  -h                  \tPrint this help message and exit."},
  {GlobalOptionIndex::Verbose, 0, "v", "verbose", validate::Verbosity,
   "  --verbose, -v <level>        \tSet verbosity level:\n"
   "                               \t  0 = errors only\n"
   "                               \t  1 = info, warnings, and errors\n"
   "                               \t  2 = debug, info, warnings, and errors"},
  {GlobalOptionIndex::ProgressMs, 0, "", "progress-ms", validate::Required,
   "  --progress-ms <ms>           \tProgress update interval in milliseconds (0 disables; default: 50)."},
  {0, 0, 0, 0, 0, 0}
};
} // namespace global_opts

[[nodiscard]] const std::unordered_map<std::string, CommandFactory>& get_command_registry() noexcept {
  static const std::unordered_map<std::string, CommandFactory> registry = {
    {"contacts",  &ContactsCommand::create},
    {"createdb",  &CreateDbCommand::create},
    {"extract",   &ExtractCommand::create},
    {"positions", &PositionsCommand::create}
  };
  return registry;
}

void print_global_help() {
  option::printUsage(std::cout, global_opts::usage);
}

} // namespace lahuta::cli
