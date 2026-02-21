/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto noop = [](const char*) {};
 *   std::unique_ptr<const char, decltype(noop)> a("besian", noop);
 *   std::unique_ptr<const char, decltype(noop)> b("sejdiu", noop);
 *   std::unique_ptr<const char, decltype(noop)> c("@gmail.com", noop);
 *   return std::string(a.get()) + b.get() + c.get();
 * }();
 *
 */

#include <pybind11/pybind11.h>

#include "pipeline/backpressure.hpp"
#include "pipeline/context.hpp"
#include "pipeline/sinks.hpp"
#include "pipeline/sources.hpp"
#include "pipeline/stage_manager.hpp"

namespace py = pybind11;
namespace lahuta::bindings {

void bind_pipeline_dynamic(py::module_ &m) {
  py::module_ md = m.def_submodule("pipeline", "Pipeline bindings");
  bind_backpressure(md);
  bind_pipeline_context(md);
  bind_sinks(md);
  bind_sources(md);
  bind_stage_manager(md);
}

} // namespace lahuta::bindings
