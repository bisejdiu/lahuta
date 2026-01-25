#ifndef LAHUTA_CLI_PIPELINE_PLAN_HPP
#define LAHUTA_CLI_PIPELINE_PLAN_HPP

#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "pipeline/dynamic/manager.hpp"

namespace lahuta::cli {

using pipeline::compute::BuildTopologyParams;
using pipeline::compute::Computation;
using pipeline::compute::Mut;
using pipeline::compute::PipelineContext;
using pipeline::compute::SystemReadParams;
using pipeline::dynamic::BackpressureConfig;
using pipeline::dynamic::IDynamicSink;
using pipeline::dynamic::ITask;
using pipeline::dynamic::StageManager;

struct PipelineReporter;

struct PipelineTask {
  std::string name;
  std::vector<std::string> deps;
  std::shared_ptr<ITask> task;
  bool thread_safe = true;
};

struct PipelineComputation {
  std::string name;
  std::vector<std::string> deps;
  std::function<std::unique_ptr<Computation<PipelineContext, Mut::ReadWrite>>()> factory;
  bool thread_safe = true;
};

struct PipelineSink {
  std::string channel;
  std::shared_ptr<IDynamicSink> sink;
  std::optional<BackpressureConfig> backpressure;
};

struct PipelinePlan {
  using SourcePtr     = StageManager::SourcePtr;
  using SourceFactory = std::function<SourcePtr()>;

  SourcePtr source;
  SourceFactory source_factory;
  bool auto_builtins = false;
  std::optional<StageManager::ReportingLevel> reporting_level;
  SystemReadParams system_params{};
  BuildTopologyParams topology_params{};
  bool override_system_params      = false;
  bool override_topology_params    = false;
  std::size_t threads              = 1;
  const PipelineReporter *reporter = nullptr;
  std::optional<std::size_t> total_items;
  std::string report_label;
  bool save_run_report = false;
  std::string run_report_prefix;
  std::string success_message;
  std::vector<PipelineTask> tasks;
  std::vector<PipelineComputation> computations;
  std::vector<PipelineSink> sinks;

  [[nodiscard]] SourcePtr resolve_source() const;
  [[nodiscard]] std::unique_ptr<StageManager> build_manager() const;

  [[nodiscard]] static PipelinePlan make_noop();
};

} // namespace lahuta::cli

#endif // LAHUTA_CLI_PIPELINE_PLAN_HPP
