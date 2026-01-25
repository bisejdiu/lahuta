#ifndef LAHUTA_CLI_COMMAND_RUNNER_HPP
#define LAHUTA_CLI_COMMAND_RUNNER_HPP

#include "specs/command_spec.hpp"

namespace lahuta::cli {

class CommandRunner {
public:
  static int run(const CommandSpec &spec, int argc, const char *const *argv);
};

} // namespace lahuta::cli

#endif // LAHUTA_CLI_COMMAND_RUNNER_HPP
