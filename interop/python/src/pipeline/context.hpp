/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto f = [](auto&&... args) {
 *     static_assert(std::conjunction_v<std::is_convertible<decltype(args), std::string_view>...>);
 *     return (std::string{} + ... + std::string(args));
 *   };
 *   return f("besian", "sejdiu", "@gmail.com");
 * }();
 *
 */

#ifndef LAHUTA_BINDINGS_PIPELINE_CONTEXT_HPP
#define LAHUTA_BINDINGS_PIPELINE_CONTEXT_HPP

#include <cstdint>
#include <string>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "numpy_utils.hpp"
#include "pipeline/data/frame.hpp"
#include "pipeline/data/model_payload.hpp"
#include "pipeline/task/context.hpp"

namespace py = pybind11;
namespace lahuta::bindings {
namespace P = lahuta::pipeline;

class PyModelPayload {
public:
  explicit PyModelPayload(std::shared_ptr<const P::ModelPayloadSlices> payload)
      : payload_(std::move(payload)) {}

  py::object sequence() const {
    if (!payload_ || !payload_->sequence) return py::none();
    return py::str(*payload_->sequence);
  }

  py::object sequence_view() const {
    if (!payload_ || !payload_->sequence_view) return py::none();
    const auto &handle = payload_->sequence_view;
    //
    // LMDB views are thread-bound. Destroying this capsule off the creating
    // thread can invalidate the txn. Copy if you need to keep data after the
    // task.
    //
    auto capsule = py::capsule(new std::shared_ptr<const P::SequenceView>(handle), [](void *p) {
      delete static_cast<std::shared_ptr<const P::SequenceView> *>(p);
    });
    auto arr     = py::array(py::dtype::of<std::uint8_t>(),
                             {static_cast<py::ssize_t>(handle->data.size())},
                             {py::ssize_t(1)},
                         reinterpret_cast<const std::uint8_t *>(handle->data.data()),
                         capsule);
    numpy::set_readonly(arr);
    return arr;
  }

  py::object plddts() const {
    if (!payload_ || !payload_->plddts || payload_->plddts->empty()) return py::none();
    const auto &cats = *payload_->plddts;
    auto arr         = py::array_t<std::uint8_t>(static_cast<py::ssize_t>(cats.size()));
    if (!cats.empty()) {
      std::memcpy(arr.mutable_data(),
                  reinterpret_cast<const std::uint8_t *>(cats.data()),
                  cats.size() * sizeof(std::uint8_t));
    }
    return arr;
  }

  py::object plddts_view() const {
    if (!payload_ || !payload_->plddts_view) return py::none();
    const auto &handle = payload_->plddts_view;
    auto capsule       = py::capsule(new std::shared_ptr<const P::PlddtView>(handle),
                               [](void *p) { delete static_cast<std::shared_ptr<const P::PlddtView> *>(p); });
    auto ptr           = reinterpret_cast<const std::uint8_t *>(handle->data.data());
    auto arr           = py::array(py::dtype::of<std::uint8_t>(),
                                   {static_cast<py::ssize_t>(handle->data.size())},
                                   {py::ssize_t(sizeof(std::uint8_t))},
                         ptr,
                         capsule);
    numpy::set_readonly(arr);
    return arr;
  }

  py::object dssp() const {
    if (!payload_ || !payload_->dssp || payload_->dssp->empty()) return py::none();
    const auto &vals = *payload_->dssp;
    auto arr         = py::array_t<std::uint8_t>(static_cast<py::ssize_t>(vals.size()));
    if (!vals.empty()) {
      std::memcpy(arr.mutable_data(),
                  reinterpret_cast<const std::uint8_t *>(vals.data()),
                  vals.size() * sizeof(std::uint8_t));
    }
    return arr;
  }

  py::object dssp_view() const {
    if (!payload_ || !payload_->dssp_view) return py::none();
    const auto &handle = payload_->dssp_view;
    auto capsule       = py::capsule(new std::shared_ptr<const P::DsspView>(handle),
                               [](void *p) { delete static_cast<std::shared_ptr<const P::DsspView> *>(p); });
    auto ptr           = reinterpret_cast<const std::uint8_t *>(handle->data.data());
    auto arr           = py::array(py::dtype::of<std::uint8_t>(),
                                   {static_cast<py::ssize_t>(handle->data.size())},
                                   {py::ssize_t(sizeof(std::uint8_t))},
                         ptr,
                         capsule);
    numpy::set_readonly(arr);
    return arr;
  }

  py::object metadata() const {
    if (!payload_ || !payload_->metadata) return py::none();
    py::dict out;
    auto emit = [&out](const char *key, const std::string &value) {
      if (value.empty()) {
        out[py::str(key)] = py::none();
      } else {
        out[py::str(key)] = py::str(value);
      }
    };
    emit("ncbi_taxonomy_id", payload_->metadata->ncbi_taxonomy_id);
    emit("organism_scientific", payload_->metadata->organism_scientific);
    return out;
  }

  py::object positions_copy() const {
    if (!payload_) return py::none();
    // Prefer decoded copy if available
    if (payload_->positions && !payload_->positions->empty()) {
      return numpy::positions_copy_f32(*payload_->positions);
    }
    // Fallback to copy from zero-copy view into a new numpy array
    if (payload_->positions_view) {
      const auto span = payload_->positions_view->data;
      auto out = py::array_t<float>({static_cast<py::ssize_t>(span.size()), static_cast<py::ssize_t>(3)});
      auto buf = out.request();
      std::memcpy(buf.ptr, span.data(), span.size() * sizeof(RDGeom::Point3Df));
      return out; // is writeable
    }
    return py::none();
  }

  py::object positions_view() const {
    if (!payload_ || !payload_->positions_view) return py::none();
    const auto &handle = payload_->positions_view;
    const auto span    = handle->data;
    const auto *ptr    = span.empty() ? nullptr : reinterpret_cast<const float *>(span.data());

    auto capsule = py::capsule(new std::shared_ptr<const P::CoordinateView>(handle), [](void *p) {
      delete static_cast<std::shared_ptr<const P::CoordinateView> *>(p);
    });

    py::array arr(py::dtype::of<float>(),
                  {static_cast<py::ssize_t>(span.size()), static_cast<py::ssize_t>(3)},
                  {static_cast<py::ssize_t>(3 * sizeof(float)), static_cast<py::ssize_t>(sizeof(float))},
                  ptr,
                  capsule);
    numpy::set_readonly(arr);
    return arr;
  }

private:
  std::shared_ptr<const P::ModelPayloadSlices> payload_;
};

class PyPipelineContext {
public:
  PyPipelineContext(P::TaskContext *ctx, std::string path) : ctx_(ctx), path_(std::move(path)) {}

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

  py::object get_system() const {
    if (!ctx_) return py::none();
    auto sys = ctx_->system();
    if (!sys) return py::none();
    return py::cast(sys);
  }

  py::object get_topology() const {
    if (!ctx_) return py::none();
    auto top = ctx_->topology();
    if (!top) return py::none();
    return py::cast(top);
  }

  //
  // Setters to enable DI builders from Python tasks
  // system and topology setters are not exposed.
  //
  void set_text(const std::string &key, const std::string &value) {
    if (!ctx_) return;
    ctx_->set_text(key, value);
  }

  void set_json(const std::string &key, py::object obj) {
    if (!ctx_) return;
    py::gil_scoped_acquire gil;
    try {
      static py::object dumps = py::module_::import("json").attr("dumps");
      std::string payload     = dumps(std::move(obj), py::arg("ensure_ascii") = false).cast<std::string>();
      ctx_->set_text(key, std::move(payload));
    } catch (const py::error_already_set &) {
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
    } catch (const py::error_already_set &) {
      // just ignore
    }
  }

  py::object get_text(const std::string &key) const {
    if (!ctx_) return py::none();
    if (auto s = ctx_->get_text(key)) return py::str(*s);
    return py::none();
  }

  py::object get_json(const std::string &key) const {
    if (!ctx_) return py::none();

    auto s = ctx_->get_text(key);
    if (!s) return py::none();

    try {
      return py::module_::import("json").attr("loads")(py::str(*s));
    } catch (const py::error_already_set &) {
      return py::none();
    }
  }

  py::object get(const std::string &key) const {
    if (!ctx_) return py::none();

    auto s = ctx_->get_text(key);
    if (!s) return py::none();

    try {
      return py::module_::import("json").attr("loads")(py::str(*s));
    } catch (const py::error_already_set &) {
      return py::str(*s);
    }
  }

  py::object get_bytes(const std::string &key) const {
    if (!ctx_) return py::none();
    if (auto s = ctx_->get_bytes(key)) return py::bytes(*s);
    return py::none();
  }

  py::object frame_metadata() const {
    if (auto meta = frame_metadata_ptr()) {
      py::dict out;
      out["session_id"]   = meta->session_id;
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

  py::object model_payload() const {
    auto payload = payload_ptr();
    if (!payload) return py::none();
    if (!cached_payload_) {
      cached_payload_ = py::cast(std::make_shared<PyModelPayload>(payload));
    }
    return cached_payload_;
  }

  py::object frame_positions_view() const {
    if (!ctx_) return py::none();
    auto conf = ctx_->conformer();
    if (!conf) return py::none();
    const auto &coords = conf->getPositions();
    if (coords.empty()) return py::none();
    auto capsule = py::capsule(new std::shared_ptr<const RDKit::Conformer>(conf), [](void *p) {
      delete static_cast<std::shared_ptr<const RDKit::Conformer> *>(p);
    });
    auto arr = numpy::make_coordinates_view_f64(coords, capsule);
    numpy::set_readonly(arr);
    return arr;
  }

private:
  const P::FrameMetadata *frame_metadata_ptr() const {
    if (!ctx_) return nullptr;
    auto meta = ctx_->frame_metadata();
    return meta ? meta.get() : nullptr;
  }

  std::shared_ptr<const P::ModelPayloadSlices> payload_ptr() const {
    if (!ctx_) return {};
    return ctx_->model_payload();
  }

  P::TaskContext *ctx_;
  std::string path_;
  mutable py::object cached_payload_;
};

inline void bind_pipeline_context(py::module_ &md) {
  // clang-format off
  constexpr const char* txn_help = R"doc( The view is backed by LMDB memory and keeps the database transaction alive while in use. Do not keep this view after your task completes. Copy the data if you need it later.)doc";
  py::class_<PyModelPayload, std::shared_ptr<PyModelPayload>>(md, "ModelPayload")
      .def_property_readonly("sequence",       &PyModelPayload::sequence)
      .def_property_readonly("sequence_view",  &PyModelPayload::sequence_view)
      .def_property_readonly("plddts",         &PyModelPayload::plddts)
      .def_property_readonly("plddts_view",    &PyModelPayload::plddts_view)
      .def_property_readonly("dssp",           &PyModelPayload::dssp)
      .def_property_readonly("dssp_view",      &PyModelPayload::dssp_view)
      .def_property_readonly("metadata",       &PyModelPayload::metadata)
      .def_property_readonly("positions",      &PyModelPayload::positions_copy)
      .def_property_readonly("positions_view", &PyModelPayload::positions_view);

  py::class_<PyPipelineContext>(md, "PipelineContext")
      .def_property_readonly("path",          &PyPipelineContext::path,          py::doc(R"doc(Return the input item path (e.g., file path).)doc"))
      .def_property_readonly("conformer_id",  &PyPipelineContext::conformer_id,  py::doc(R"doc(Return the conformer/frame identifier for the current item.)doc"))
      .def_property_readonly("session_id",    &PyPipelineContext::session_id,    py::doc(R"doc(Return the session identifier associated with this item.)doc"))
      .def_property_readonly("timestamp_ps",  &PyPipelineContext::timestamp_ps,  py::doc(R"doc(Optional simulation timestamp in picoseconds.)doc"))
      .def_property_readonly("model_payload", &PyPipelineContext::model_payload, py::doc(R"doc(Structured view over lazily loaded LMDB payload slices.)doc"))
      .def_property_readonly("frame_positions_view", &PyPipelineContext::frame_positions_view, py::doc(R"doc(Active frame coordinates as a read-only NumPy view with shape (n_atoms, 3).)doc"))

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
