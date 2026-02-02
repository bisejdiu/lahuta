#ifndef LAHUTA_CLI_TASKS_COMPACTION_RG_SUMMARY_SINK_HPP
#define LAHUTA_CLI_TASKS_COMPACTION_RG_SUMMARY_SINK_HPP

#include <filesystem>
#include <memory>

#include "pipeline/io/sink_iface.hpp"
#include "tasks/compaction_rg_task.hpp"

namespace lahuta::cli::compaction_rg {
namespace P = lahuta::pipeline;

[[nodiscard]] std::shared_ptr<P::IDynamicSink>
make_compaction_rg_summary_sink(std::filesystem::path output_path,
                                std::shared_ptr<CompactionRgCounters> counters);

} // namespace lahuta::cli::compaction_rg

#endif // LAHUTA_CLI_TASKS_COMPACTION_RG_SUMMARY_SINK_HPP
