#ifndef LAHUTA_CLI_TASKS_COMPACTION_RG_TASK_HPP
#define LAHUTA_CLI_TASKS_COMPACTION_RG_TASK_HPP

#include <memory>
#include <string_view>

#include "analysis/compaction/rg_utils.hpp"
#include "pipeline/task/task.hpp"
#include "tasks/processing_counters.hpp"

namespace lahuta::cli::compaction_rg {
namespace P = lahuta::pipeline;

using CompactionRgCounters = ProcessingCounters<analysis::ThresholdCount>;

struct CompactionRgConfig {
  double min_high_fraction = 0.80;
  std::shared_ptr<CompactionRgCounters> counters;
};

constexpr std::string_view OutputChannel = "per_protein_rg";

[[nodiscard]] std::shared_ptr<P::ITask>
make_compaction_rg_task(std::shared_ptr<const CompactionRgConfig> config);

} // namespace lahuta::cli::compaction_rg

#endif // LAHUTA_CLI_TASKS_COMPACTION_RG_TASK_HPP
