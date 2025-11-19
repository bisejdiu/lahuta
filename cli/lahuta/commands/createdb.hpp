#ifndef LAHUTA_CLI_CREATEDB_HPP
#define LAHUTA_CLI_CREATEDB_HPP

#include <memory>

#include "commands/command.hpp"

namespace lahuta::cli {

namespace createdb_opts {
enum CreateDbOptionIndex : unsigned {
  Unknown,
  Help,

  // Source options
  SourceDirectory,
  SourceVector,
  SourceFileList,
  Extension,
  Recursive,

  // Database options
  DatabasePath,
  BatchSize,
  MaxSize,

  // Runtime options
  Threads
};
} // namespace createdb_opts

class CreateDbCommand final : public CliCommand {
public:
  [[nodiscard]] static std::unique_ptr<CliCommand> create();
  int run(int argc, char* argv[]) override;

private:
  CreateDbCommand() = default;
};

} // namespace lahuta::cli

#endif // LAHUTA_CLI_CREATEDB_HPP
