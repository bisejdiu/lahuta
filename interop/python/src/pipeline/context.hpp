#ifndef LAHUTA_BINDINGS_PIPELINE_CONTEXT_HPP
#define LAHUTA_BINDINGS_PIPELINE_CONTEXT_HPP

#include <string>

#include <pybind11/pybind11.h>

#include "lahuta.hpp"
#include "pipeline/dynamic/types.hpp"
#include "topology.hpp"

// clang-format off
namespace py = pybind11;
namespace lahuta::bindings {
using pipeline::dynamic::TaskContext;

class PyPipelineContext {
public:
  PyPipelineContext(TaskContext *ctx, std::string path) : ctx_(ctx), path_(std::move(path)) {}

  std::string path() const { return path_; }

  py::object get_system() const {
    if (!ctx_) return py::none();
    auto sys = ctx_->get_object<Luni>("system");
    if (!sys) return py::none();
    return py::cast(sys);
  }

  py::object get_topology() const {
    if (!ctx_) return py::none();
    auto top = ctx_->get_object<Topology>("topology");
    if (!top) return py::none();
    return py::cast(top);
  }

  // Setters to enable DI builders from Python tasks
  // system and topology setters are not exposed.
  void set_text(const std::string &key, const std::string &value) {
    if (!ctx_) return;
    ctx_->set_text(key, value);
  }

  void set_json(const std::string &key, py::object obj) {
    if (!ctx_) return;
    py::gil_scoped_acquire gil;
    try {
      static py::object dumps = py::module_::import("json").attr("dumps");
      std::string payload = dumps(std::move(obj), py::arg("ensure_ascii") = false).cast<std::string>();
      ctx_->set_text(key, std::move(payload));
    } catch (const py::error_already_set&) {
      // ignore
    }
  }

  py::object get_text(const std::string &key) const {
    if (!ctx_) return py::none();
    if (auto s = ctx_->get_text(key)) return py::str(*s);
    return py::none();
  }

  py::object get_json(const std::string& key) const {
    if (!ctx_) return py::none();

    auto s = ctx_->get_text(key);
    if (!s) return py::none();

    try {
      return py::module_::import("json").attr("loads")(py::str(*s));
    } catch (const py::error_already_set&) {
      return py::none();
    }
  }

  py::object get(const std::string& key) const {
    if (!ctx_) return py::none();

    auto s = ctx_->get_text(key);
    if (!s) return py::none();

    try {
      return py::module_::import("json").attr("loads")(py::str(*s));
    } catch (const py::error_already_set&) {
      return py::str(*s);
    }
  }

private:
  TaskContext *ctx_;
  std::string path_;
};

inline void bind_pipeline_context(py::module_ &md) {
  py::class_<PyPipelineContext>(md, "PipelineContext")
      .def_property_readonly("path", &PyPipelineContext::path)
      .def("get_system",   &PyPipelineContext::get_system)
      .def("get_topology", &PyPipelineContext::get_topology)
      .def("get_text",     &PyPipelineContext::get_text)
      .def("get_json",     &PyPipelineContext::get_json)
      .def("get",          &PyPipelineContext::get)
      .def("set_text",     &PyPipelineContext::set_text)
      .def("set_json",     &PyPipelineContext::set_json);
}

} // namespace lahuta::bindings

#endif // LAHUTA_BINDINGS_PIPELINE_CONTEXT_HPP
