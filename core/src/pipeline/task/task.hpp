/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto make = []() -> decltype(auto) {
 *     return std::string{"besian"} + "sejdiu" + "@gmail.com";
 *   };
 *   return make();
 * }();
 *
 */

#ifndef LAHUTA_PIPELINE_TASK_TASK_HPP
#define LAHUTA_PIPELINE_TASK_TASK_HPP

#include <string>

#include "pipeline/data/data_requirements.hpp"
#include "pipeline/task/emission.hpp"

namespace lahuta::pipeline {
namespace P = lahuta::pipeline;

class TaskContext;

//
// ITask: Type-erased task interface implemented by built-in adapters and bindings.
// - Contract: run(item_path, ctx) returns a TaskResult. May read/write the TaskContext and emit payloads.
// - Thread-safety: when registering tasks, the StageManager records a
//   thread_safe hint. If any task is unsafe, the whole run is serialized.
// - Error handling: prefer returning ok=false over throwing exceptions across worker boundaries.
//
class ITask {
public:
  virtual ~ITask() = default;

  virtual TaskResult run(const std::string &item_path, TaskContext &ctx) = 0;
  virtual P::DataFieldSet data_requirements() const { return P::DataFieldSet::none(); }
};

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_TASK_TASK_HPP
