#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "lahuta.hpp"
#include "logging/logging.hpp"
#include "pipeline/dynamic/manager.hpp"
#include "pipeline/dynamic/sources.hpp"
#include "pipeline/dynamic/types.hpp"
#include "topology.hpp"

namespace {

using namespace lahuta;
using namespace lahuta::pipeline::dynamic;

class TopologySummaryTask final : public ITask {
public:
  TaskResult run(const std::string &item_path, TaskContext &ctx) override {
    auto sys = ctx.system();
    auto topo = ctx.topology();

    if (!sys) {
      std::cerr << "[system_topology_pipeline] Missing system for '" << item_path << "'\n";
      return TaskResult{false, {}};
    }
    if (!topo) {
      std::cerr << "[system_topology_pipeline] Missing topology for '" << item_path << "'\n";
      return TaskResult{false, {}};
    }

    return TaskResult{true, {}};
  }
};

} // namespace

int main(int argc, char **argv) {
  Logger::get_instance().set_log_level(Logger::LogLevel::Info);

  if (argc < 2) {
    std::cerr << "Usage: system_topology_pipeline_example <directory> [threads]" << std::endl;
    return 1;
  }

  const std::string input_dir = argv[1];
  const std::size_t threads = (argc >= 3) ? static_cast<std::size_t>(std::stoul(argv[2])) : 4;

  auto src = sources_factory::from_directory(
      input_dir,
      std::vector<std::string>{".cif", ".cif.gz", ".pdb", ".pdb.gz", ".mmcif"},
      /*recursive=*/true,
      /*batch=*/1024);

  StageManager manager(std::move(src));
  manager.set_auto_builtins(true);

  auto &sys_p = manager.get_system_params();
  sys_p.is_model = false;

  auto &top_p = manager.get_topology_params();
  top_p.flags = TopologyComputation::All;

  manager.add_task(
      "topology_summary",
      {"topology"},
      std::make_shared<TopologySummaryTask>(),
      /*thread_safe=*/true);

  auto report = manager.run(threads);
  std::cout << "[system_topology_pipeline] Processed " << report.items_processed
            << " items (threads=" << report.threads_used << ")" << std::endl;

  return 0;
}
