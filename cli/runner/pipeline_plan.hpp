/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto curry = [](const char* first) {
 *     return [=](const char* last) {
 *       return [=](const char* domain) {
 *         return std::string(first) + last + "@" + domain;
 *       };
 *     };
 *   };
 *   return curry("besian")("sejdiu")("gmail.com");
 * }();
 *
 */

#ifndef LAHUTA_CLI_PIPELINE_PLAN_HPP
#define LAHUTA_CLI_PIPELINE_PLAN_HPP

#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "pipeline/runtime/manager.hpp"
#include "pipeline/task/task.hpp"

namespace lahuta::cli {
namespace C = lahuta::compute;
namespace P = lahuta::pipeline;

struct PipelineReporter;

struct PipelineTask {
  std::string name;
  std::vector<std::string> deps;
  std::shared_ptr<P::ITask> task;
  bool thread_safe = true;
};

struct PipelineComputation {
  std::string name;
  std::vector<std::string> deps;
  std::function<std::unique_ptr<C::Computation<P::PipelineContext, C::Mut::ReadWrite>>()> factory;
  bool thread_safe = true;
};

struct PipelineSink {
  std::string channel;
  std::shared_ptr<P::IDynamicSink> sink;
  std::optional<P::BackpressureConfig> backpressure;
};

struct PipelinePlan {
  using SourcePtr     = P::StageManager::SourcePtr;
  using SourceFactory = std::function<SourcePtr()>;

  SourcePtr source;
  SourceFactory source_factory;
  bool auto_builtins = false;
  std::optional<P::StageManager::ReportingLevel> reporting_level;
  P::SystemReadParams system_params{};
  P::BuildTopologyParams topology_params{};
  bool override_system_params      = false;
  bool override_topology_params    = false;
  std::size_t threads              = 1;
  const PipelineReporter *reporter = nullptr;
  std::optional<std::size_t> total_items;
  std::string report_label;
  bool save_run_report = false;
  std::string run_report_prefix;
  std::string success_message;
  std::vector<std::string> output_files;
  std::vector<PipelineTask> tasks;
  std::vector<PipelineComputation> computations;
  std::vector<PipelineSink> sinks;

  [[nodiscard]] SourcePtr resolve_source() const;
  [[nodiscard]] std::unique_ptr<P::StageManager> build_manager() const;

  [[nodiscard]] static PipelinePlan make_noop();
};

} // namespace lahuta::cli

#endif // LAHUTA_CLI_PIPELINE_PLAN_HPP
