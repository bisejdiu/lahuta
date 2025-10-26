#ifndef LAHUTA_BINDINGS_PIPELINE_CONTEXT_HPP
#define LAHUTA_BINDINGS_PIPELINE_CONTEXT_HPP

#include <cstdint>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "lahuta.hpp"
#include "models/dssp.hpp"
#include "models/metadata.hpp"
#include "models/plddt.hpp"
#include "pipeline/dynamic/keys.hpp"
#include "pipeline/dynamic/types.hpp"
#include "pipeline/frame.hpp"
#include "topology.hpp"

// clang-format off
namespace py = pybind11;
namespace lahuta::bindings {
using pipeline::dynamic::TaskContext;

class PyPipelineContext {
public:
  PyPipelineContext(TaskContext *ctx, std::string path) : ctx_(ctx), path_(std::move(path)) {}

  std::string path() const { return path_; }

  std::uint64_t conformer_id() const {
    if (auto meta = frame_metadata_ptr()) {
      return meta->conformer_id;
    }
    return 0;
  }

  py::object session_id() const {
    if (auto meta = frame_metadata_ptr()) {
      if (!meta->session_id.empty()) {
        return py::str(meta->session_id);
      }
    }
    return py::none();
  }

  py::object timestamp_ps() const {
    if (auto meta = frame_metadata_ptr()) {
      if (meta->timestamp_ps.has_value()) {
        return py::float_(*meta->timestamp_ps);
      }
    }
    return py::none();
  }

  py::object taxonomy_metadata() const {
    if (!ctx_) return py::none();
    auto meta = ctx_->get_object<ModelMetadata>(pipeline::CTX_MODEL_METADATA_KEY);
    if (!meta) return py::none();
    py::dict out;
    if (meta->ncbi_taxonomy_id.empty()) {
      out["ncbi_taxonomy_id"] = py::none();
    } else {
      out["ncbi_taxonomy_id"] = py::str(meta->ncbi_taxonomy_id);
    }
    if (meta->organism_scientific.empty()) {
      out["organism_scientific"] = py::none();
    } else {
      out["organism_scientific"] = py::str(meta->organism_scientific);
    }
    return out;
  }

  py::object get_system() const {
    if (!ctx_) return py::none();
    auto sys = ctx_->get_object<const Luni>("system");
    if (!sys) return py::none();
    return py::cast(sys);
  }

  py::object get_topology() const {
    if (!ctx_) return py::none();
    auto top = ctx_->get_object<const Topology>("topology");
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

  void set_bytes(const std::string &key, py::object obj) {
    if (!ctx_) return;
    py::gil_scoped_acquire gil;
    try {
      // Accept bytes-like (bytes, bytearray, memoryview), py::bytes will coerce
      py::bytes b = py::bytes(obj);
      ctx_->set_bytes(key, b.cast<std::string>());
    } catch (const py::error_already_set&) {
      // invalid types are ignored
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

  py::object get_bytes(const std::string &key) const {
    if (!ctx_) return py::none();
    if (auto s = ctx_->get_bytes(key)) return py::bytes(*s);
    return py::none();
  }

  py::object plddt_categories() const {
    if (!ctx_) return py::none();
    auto cats = ctx_->get_object<const std::vector<pLDDTCategory>>(pipeline::CTX_PLDDT_KEY);
    if (!cats) return py::none();
    py::list out(cats->size());
    for (size_t i = 0; i < cats->size(); ++i) {
      out[i] = py::int_(static_cast<std::uint8_t>((*cats)[i]));
    }
    return out;
  }

  py::object dssp_assignments() const {
    if (!ctx_) return py::none();
    auto dssp = ctx_->get_object<const std::vector<DSSPAssignment>>(pipeline::CTX_DSSP_KEY);
    if (!dssp) return py::none();
    py::list out(dssp->size());
    for (size_t i = 0; i < dssp->size(); ++i) {
      out[i] = py::int_(static_cast<std::uint8_t>((*dssp)[i]));
    }
    return out;
  }

  py::object frame_metadata() const {
    if (auto meta = frame_metadata_ptr()) {
      py::dict out;
      out["session_id"] = meta->session_id;
      out["conformer_id"] = py::int_(meta->conformer_id);
      if (meta->timestamp_ps.has_value()) {
        out["timestamp_ps"] = py::float_(*meta->timestamp_ps);
      } else {
        out["timestamp_ps"] = py::none();
      }
      return out;
    }
    return py::none();
  }

private:
  const FrameMetadata* frame_metadata_ptr() const {
    if (!ctx_) return nullptr;
    auto meta = ctx_->get_object<FrameMetadata>("lahuta.frame");
    return meta ? meta.get() : nullptr;
  }

  TaskContext *ctx_;
  std::string path_;
};

inline void bind_pipeline_context(py::module_ &md) {
  py::class_<PyPipelineContext>(md, "PipelineContext")
      .def_property_readonly("path",         &PyPipelineContext::path,         py::doc(R"doc(Return the input item path (e.g., file path).)doc"))
      .def_property_readonly("conformer_id", &PyPipelineContext::conformer_id, py::doc(R"doc(Return the conformer/frame identifier for the current item.)doc"))
      .def_property_readonly("session_id",   &PyPipelineContext::session_id,   py::doc(R"doc(Return the session identifier associated with this item.)doc"))
      .def_property_readonly("timestamp_ps", &PyPipelineContext::timestamp_ps, py::doc(R"doc(Optional simulation timestamp in picoseconds.)doc"))
      .def_property_readonly("taxonomy",     &PyPipelineContext::taxonomy_metadata, py::doc(R"doc(Model taxonomy metadata for AlphaFold-like inputs (dict with ncbi_taxonomy_id / organism_scientific).)doc"))
      .def_property_readonly("plddt",        &PyPipelineContext::plddt_categories, py::doc(R"doc(Residue-level pLDDT categories (list of ints, 0=VeryHigh ... 3=VeryLow).)doc"))
      .def_property_readonly("dssp",         &PyPipelineContext::dssp_assignments, py::doc(R"doc(Residue-level DSSP assignments (list of ints: 0=Coil, 1=Alpha helix, 2=3-10 helix, 3=Pi helix, 4=PolyPro helix, 5=Strand, 6=Turn, 7=Bend).)doc"))

      .def("get",          &PyPipelineContext::get)
      .def("get_system",   &PyPipelineContext::get_system)
      .def("get_topology", &PyPipelineContext::get_topology)
      .def("get_text",     &PyPipelineContext::get_text)
      .def("get_json",     &PyPipelineContext::get_json)
      .def("get_bytes",    &PyPipelineContext::get_bytes)
      .def("set_text",     &PyPipelineContext::set_text)
      .def("set_json",     &PyPipelineContext::set_json)
      .def("set_bytes",    &PyPipelineContext::set_bytes)

      .def("frame_metadata", &PyPipelineContext::frame_metadata, py::doc(R"doc(Current frame metadata.)doc"));
}

} // namespace lahuta::bindings

#endif // LAHUTA_BINDINGS_PIPELINE_CONTEXT_HPP
