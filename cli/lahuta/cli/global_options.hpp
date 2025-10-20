#ifndef LAHUTA_CLI_GLOBAL_OPTIONS_HPP
#define LAHUTA_CLI_GLOBAL_OPTIONS_HPP

#include <memory>
#include <string>
#include <unordered_map>

#include "commands/command.hpp"
#include "logging.hpp"

namespace lahuta::cli {

//
// Ideally, we'd want to use scoped enums here, but then we'd have to explicitly cast every time
// we'd want to use them, which is a lot of places. We we use plain enums wrapped in a namespace instead.
//
namespace global_opts {
enum GlobalOptionIndex : unsigned {
  Unknown,
  Help,
  Verbose
};
} // namespace global_opts

using CommandFactory = std::unique_ptr<CliCommand>(*)();

// Get the command registry mapping subcommand names to factory functions
[[nodiscard]] const std::unordered_map<std::string, CommandFactory>& get_command_registry() noexcept;

// Parse global options and return the subcommand name and remaining arguments. Returns subcommand name if found, or empty string on error/help
[[nodiscard]] std::string parse_global_options(int argc, char* argv[], lahuta::Logger::LogLevel& log_level, int& sub_argc, char**& sub_argv);

void print_global_help();

} // namespace lahuta::cli

#endif // LAHUTA_CLI_GLOBAL_OPTIONS_HPP
