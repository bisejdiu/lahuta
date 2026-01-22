#ifndef LAHUTA_CLI_EXTRACT_HPP
#define LAHUTA_CLI_EXTRACT_HPP

#include <memory>

#include "commands/command.hpp"

namespace lahuta::cli {

namespace extract_opts {
enum ExtractOptionIndex : unsigned {
  Unknown,
  Help,
  SourceDatabase,
  SourceDirectory,
  SourceVector,
  SourceFileList,
  Extension,
  Recursive,
  Output,
  Reporter,
  SaveRunReport,
  Threads,
  BatchSize,
  WriterThreads,
  IsAf2Model
};
} // namespace extract_opts

class ExtractCommand final : public CliCommand {
public:
  [[nodiscard]] static std::unique_ptr<CliCommand> create();
  int run(int argc, char *argv[]) override;

private:
  ExtractCommand() = default;
};

} // namespace lahuta::cli

#endif // LAHUTA_CLI_EXTRACT_HPP
