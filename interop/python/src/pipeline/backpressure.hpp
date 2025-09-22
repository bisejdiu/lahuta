#ifndef LAHUTA_BINDINGS_PIPELINE_BACKPRESSURE_HPP
#define LAHUTA_BINDINGS_PIPELINE_BACKPRESSURE_HPP

#include <pybind11/chrono.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pipeline/dynamic/backpressure.hpp"

namespace py = pybind11;
namespace lahuta::bindings {
using namespace lahuta::pipeline::dynamic;

inline void bind_backpressure(py::module_& md) {
  py::enum_<OnFull>(md, "OnFull", "Policy when a sink queue is full")
      .value("Block", OnFull::Block)
      .value("DropLatest", OnFull::DropLatest)
      .value("DropOldest", OnFull::DropOldest)
      .export_values();

  py::class_<BackpressureConfig>(md, "BackpressureConfig", "Sink backpressure configuration")
      .def(py::init<>())
      .def_readwrite("max_queue_msgs",  &BackpressureConfig::max_queue_msgs)
      .def_readwrite("max_queue_bytes", &BackpressureConfig::max_queue_bytes)
      .def_readwrite("max_batch_msgs",  &BackpressureConfig::max_batch_msgs)
      .def_readwrite("max_batch_bytes", &BackpressureConfig::max_batch_bytes)
      .def_readwrite("offer_timeout",   &BackpressureConfig::offer_timeout,
          "Producer wait slice between retries when the sink queue is full. It's nnot a hard timeout.")
      .def_readwrite("on_full",         &BackpressureConfig::on_full)
      .def_readwrite("required",        &BackpressureConfig::required)
      .def("validate", [](const BackpressureConfig& cfg) { validate_config(cfg); });

  md.def(
      "get_default_backpressure_config",
      []() { return get_default_backpressure_config(); },
      "Return a copy of the default sink backpressure configuration."
  );

  md.def(
      "set_default_backpressure_config",
      [](const BackpressureConfig& cfg) { set_default_backpressure_config(cfg); },
      "Set the default sink backpressure configuration used for new sinks."
  );

  md.def(
      "set_default_max_queue_bytes",
      [](std::size_t bytes) { set_default_max_queue_bytes(bytes); },
      py::arg("bytes"),
      "Update only the byte capacity of the default sink queue (clamps batch bytes accordingly)."
  );
}

} // namespace lahuta::bindings

#endif // LAHUTA_BINDINGS_PIPELINE_BACKPRESSURE_HPP
