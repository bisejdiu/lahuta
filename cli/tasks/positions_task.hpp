#ifndef LAHUTA_CLI_TASKS_POSITIONS_TASK_HPP
#define LAHUTA_CLI_TASKS_POSITIONS_TASK_HPP

#include <filesystem>
#include <memory>

#include "pipeline/task/task.hpp"

namespace lahuta::cli::positions {
namespace P = lahuta::pipeline;

[[nodiscard]] std::shared_ptr<P::ITask> make_positions_task(std::filesystem::path output_dir, int tree_depth);

} // namespace lahuta::cli::positions

#endif // LAHUTA_CLI_TASKS_POSITIONS_TASK_HPP
