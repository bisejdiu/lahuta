#ifndef LAHUTA_CLI_QUALITY_METRICS_HPP
#define LAHUTA_CLI_QUALITY_METRICS_HPP

#include <memory>

#include "commands/command.hpp"

namespace lahuta::cli {

namespace quality_metrics_opts {
enum QualityMetricsOptionIndex : unsigned {
  Unknown,
  Help,
  PlddtGroup,
  DsspGroup,
  SegmentGroup,
  SegmentMin,
  NoOverlap,
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
} // namespace quality_metrics_opts

class QualityMetricsCommand final : public CliCommand {
public:
  [[nodiscard]] static std::unique_ptr<CliCommand> create();
  int run(int argc, char *argv[]) override;

private:
  QualityMetricsCommand() = default;
};

} // namespace lahuta::cli

#endif // LAHUTA_CLI_QUALITY_METRICS_HPP
