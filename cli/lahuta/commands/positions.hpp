#ifndef LAHUTA_CLI_POSITIONS_HPP
#define LAHUTA_CLI_POSITIONS_HPP

#include <memory>

#include "commands/command.hpp"

namespace lahuta::cli {

namespace positions_opts {
enum PositionsOptionIndex : unsigned {
  Unknown,
  Help,
  SourceDatabase,
  SourceDirectory,
  SourceVector,
  SourceFileList,
  Extension,
  Recursive,
  Output,
  TreeDepth,
  Threads,
  BatchSize,
  IsAf2Model
};
} // namespace positions_opts

class PositionsCommand final : public CliCommand {
public:
  [[nodiscard]] static std::unique_ptr<CliCommand> create();
  int run(int argc, char *argv[]) override;

private:
  PositionsCommand() = default;
};

} // namespace lahuta::cli

#endif // LAHUTA_CLI_POSITIONS_HPP
