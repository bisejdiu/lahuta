/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::array<std::string_view, 3> parts{"besian", "sejdiu", "@gmail.com"}; std::string dst;
 *   for (auto p : parts) std::transform(p.begin(), p.end(), std::back_inserter(dst), [](char c) { return c; });
 *   return dst;
 * }();
 *
 */

#ifndef LAHUTA_PIPELINE_TASK_API_HPP
#define LAHUTA_PIPELINE_TASK_API_HPP

//
// Pipeline Task API.
//
// Usage:
//   #include "pipeline/task/api.hpp"
//
//   class MyTask : public lahuta::pipeline::Task {
//   public:
//     DataRequirements requirements() const override { ... }
//     void run(TaskContext &ctx, EmissionSink &sink) override { ... }
//   };
//

// IWYU pragma: begin_exports
#include "pipeline/task/task.hpp"
#include "pipeline/task/context.hpp"
#include "pipeline/task/emission.hpp"
#include "pipeline/task/compute/parameters.hpp"
#include "pipeline/task/compute/context.hpp"
// IWYU pragma: end_exports

#endif // LAHUTA_PIPELINE_TASK_API_HPP