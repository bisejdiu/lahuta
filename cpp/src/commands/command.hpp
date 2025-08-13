#ifndef LAHUTA_CLI_COMMAND_HPP
#define LAHUTA_CLI_COMMAND_HPP

namespace lahuta::cli {

class CliCommand {
public:
  virtual ~CliCommand() = default;
  virtual int run(int argc, char* argv[]) = 0;
};

} // namespace lahuta::cli

#endif // LAHUTA_CLI_COMMAND_HPP
