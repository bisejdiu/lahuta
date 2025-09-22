#include <pybind11/pybind11.h>

#include "pipeline/backpressure.hpp"
#include "pipeline/context.hpp"
#include "pipeline/sinks.hpp"
#include "pipeline/stage_manager.hpp"

namespace py = pybind11;
namespace lahuta::bindings {

void bind_pipeline_dynamic(py::module_ &m) {
  py::module_ md = m.def_submodule("pipeline", "Pipeline bindings");
  bind_backpressure(md);
  bind_pipeline_context(md);
  bind_sinks(md);
  bind_stage_manager(md);
}

} // namespace lahuta::bindings
