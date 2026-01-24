#ifndef LAHUTA_CLI_GLOBAL_OPTIONS_HPP
#define LAHUTA_CLI_GLOBAL_OPTIONS_HPP

#include <memory>
#include <string>
#include <unordered_map>

#include "commands/command.hpp"
#include "logging/logging.hpp"

namespace lahuta::cli {

//
// Ideally, we'd want to use scoped enums here, but then we'd have to explicitly cast every time
// we'd want to use them, which is a lot of places. We we use plain enums wrapped in a namespace instead.
//
namespace global_opts {
enum GlobalOptionIndex : unsigned {
  Unknown,
  Help,
  Verbose,
  ProgressMs,
  ProgressNoColor
};
} // namespace global_opts

using CommandFactory = std::unique_ptr<CliCommand>(*)();

// Get the command registry mapping subcommand names to factory functions
[[nodiscard]] const std::unordered_map<std::string, CommandFactory>& get_command_registry() noexcept;

void print_global_help();

} // namespace lahuta::cli

#endif // LAHUTA_CLI_GLOBAL_OPTIONS_HPP
