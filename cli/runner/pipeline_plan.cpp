#include <stdexcept>

#include "runner/pipeline_plan.hpp"
#include "pipeline/dynamic/sources.hpp"

namespace lahuta::cli {
namespace dyn = pipeline::dynamic;

PipelinePlan::SourcePtr PipelinePlan::resolve_source() const {
  if (source_factory) {
    return source_factory();
  }
  return source;
}

std::unique_ptr<dyn::StageManager> PipelinePlan::build_manager() const {
  auto resolved = resolve_source();
  if (!resolved) {
    throw std::runtime_error("PipelinePlan requires a source descriptor");
  }

  auto manager = std::make_unique<dyn::StageManager>(std::move(resolved));
  manager->set_auto_builtins(auto_builtins);

  if (reporting_level.has_value()) {
    manager->set_reporting_level(*reporting_level);
  }

  if (override_system_params) {
    manager->get_system_params() = system_params;
  }
  if (override_topology_params) {
    manager->get_topology_params() = topology_params;
  }

  for (const auto &task : tasks) {
    if (!task.task) {
      throw std::runtime_error("PipelinePlan task '" + task.name + "' has null implementation");
    }
    manager->add_task(task.name, task.deps, task.task, task.thread_safe);
  }

  for (const auto &comp : computations) {
    if (!comp.factory) {
      throw std::runtime_error("PipelinePlan computation '" + comp.name + "' has null factory");
    }
    manager->add_computation(comp.name, comp.deps, comp.factory, comp.thread_safe);
  }

  for (const auto &sink : sinks) {
    if (!sink.sink) {
      throw std::runtime_error("PipelinePlan sink for channel '" + sink.channel + "' is null");
    }
    manager->connect_sink(sink.channel, sink.sink, sink.backpressure);
  }

  return manager;
}

PipelinePlan PipelinePlan::make_noop() {
  PipelinePlan plan;
  plan.source          = std::shared_ptr<sources::IDescriptor>(dyn::sources_factory::from_vector({}));
  plan.reporting_level = dyn::StageManager::ReportingLevel::Basic;
  plan.threads         = 1;
  return plan;
}

} // namespace lahuta::cli
