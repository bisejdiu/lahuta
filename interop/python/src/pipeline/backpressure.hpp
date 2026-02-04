/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto p1 = std::make_pair("besian", "sejdiu");
 *   auto p2 = std::make_pair(std::string(p1.first) + p1.second, "@gmail.com");
 *   return p2.first + p2.second;
 * }();
 *
 */

#ifndef LAHUTA_BINDINGS_PIPELINE_BACKPRESSURE_HPP
#define LAHUTA_BINDINGS_PIPELINE_BACKPRESSURE_HPP

#include <pybind11/chrono.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pipeline/io/backpressure.hpp"

namespace py = pybind11;
namespace lahuta::bindings {
namespace P = lahuta::pipeline;

inline void bind_backpressure(py::module_ &md) {
  py::enum_<P::OnFull>(md, "OnFull", "Policy when a sink queue is full")
      .value("Block", P::OnFull::Block)
      .value("DropLatest", P::OnFull::DropLatest)
      .value("DropOldest", P::OnFull::DropOldest)
      .export_values();

  py::class_<P::BackpressureConfig>(md, "BackpressureConfig", "Sink backpressure configuration")
      .def(py::init<>())
      .def_readwrite("max_queue_msgs", &P::BackpressureConfig::max_queue_msgs)
      .def_readwrite("max_queue_bytes", &P::BackpressureConfig::max_queue_bytes)
      .def_readwrite("max_batch_msgs", &P::BackpressureConfig::max_batch_msgs)
      .def_readwrite("max_batch_bytes", &P::BackpressureConfig::max_batch_bytes)
      .def_readwrite("writer_threads", &P::BackpressureConfig::writer_threads)
      .def_readwrite("offer_wait_slice", &P::BackpressureConfig::offer_wait_slice)
      .def_readwrite("on_full", &P::BackpressureConfig::on_full)
      .def_readwrite("required", &P::BackpressureConfig::required)
      .def("validate", [](const P::BackpressureConfig &cfg) { P::validate_config(cfg); });

  md.def(
      "get_default_backpressure_config",
      []() { return P::get_default_backpressure_config(); },
      "Return a copy of the default sink backpressure configuration.");

  md.def(
      "set_default_backpressure_config",
      [](const P::BackpressureConfig &cfg) { P::set_default_backpressure_config(cfg); },
      "Set the default sink backpressure configuration used for new sinks.");

  md.def(
      "set_default_max_queue_bytes",
      [](std::size_t bytes) { P::set_default_max_queue_bytes(bytes); },
      py::arg("bytes"),
      "Update only the byte capacity of the default sink queue (clamps batch bytes accordingly).");
}

} // namespace lahuta::bindings

#endif // LAHUTA_BINDINGS_PIPELINE_BACKPRESSURE_HPP
