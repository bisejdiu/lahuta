#ifndef LAHUTA_PIPELINE_TASKS_TOPOLOGY_TASK_HPP
#define LAHUTA_PIPELINE_TASKS_TOPOLOGY_TASK_HPP

#include "lahuta.hpp"
#include "logging.hpp"
#include <string>
#include <string_view>

// clang-format off
namespace lahuta::tasks {

class TopologyTask {
public:
  struct TestTopData {
    size_t   num_atoms;
    size_t   num_bonds;
  };

  using input_type = std::string_view;
  struct result_type {
    bool success;
    std::string file_path;
    TestTopData data;
  };

  explicit TopologyTask(TopologyComputation flags = TopologyComputation::All)
    : flags_(flags)
  {}

  result_type operator()(std::string_view file_path) const {
    result_type result;
    result.file_path = std::string(file_path);
    try {
      Luni luni(result.file_path);
      luni.enable_topology_only(flags_);
      result.success   = luni.build_topology();
      result.data.num_atoms = luni.get_molecule().getNumAtoms();
      result.data.num_bonds = luni.get_molecule().getNumBonds();
    } catch (const std::exception &e) {
      result.success = false;
      Logger::get_logger()->error(
        "Error processing file {}: {}",
        file_path, e.what()
      );
    }
    return result;
  }

private:
  TopologyComputation flags_;
};

} // namespace lahuta::tasks

#endif // LAHUTA_PIPELINE_TASKS_TOPOLOGY_TASK_HPP
