#ifndef LAHUTA_BINDINGS_DYNAMIC_SINKS_HPP
#define LAHUTA_BINDINGS_DYNAMIC_SINKS_HPP

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "io/sinks/memory.hpp"
#include "io/sinks/ndjson.hpp"
#include "io/sinks/sharded_ndjson.hpp"
#include "io/sinks/lmdb.hpp"
#include "db/db.hpp"

namespace py = pybind11;
namespace lahuta::bindings {
using namespace lahuta::pipeline::dynamic;

// clang-format off
inline void bind_sinks(py::module_& md) {
  py::class_<IDynamicSink, std::shared_ptr<IDynamicSink>> _(md, "Sink");

  py::class_<MemorySink, IDynamicSink, std::shared_ptr<MemorySink>>(md, "MemorySink")
    .def(py::init<>())
    .def("result", &MemorySink::result, R"doc(Return collected payloads.)doc")
    .def("clear",  &MemorySink::clear,  R"doc(Clear collected payloads.)doc");

  py::class_<NdjsonFileSink, IDynamicSink, std::shared_ptr<NdjsonFileSink>>(md, "NdjsonSink")
    .def(py::init<std::string>(), py::arg("path"))
    .def("file", &NdjsonFileSink::file);

  py::class_<ShardedNdjsonSink, IDynamicSink, std::shared_ptr<ShardedNdjsonSink>>(md, "ShardedNdjsonSink")
    .def(py::init<std::string, std::size_t>(), py::arg("out_dir"), py::arg("shard_size"))
    .def(py::init<std::string, std::size_t, std::size_t>(), py::arg("out_dir"), py::arg("shard_size"), py::arg("max_shard_bytes"))
    .def("files", &ShardedNdjsonSink::files);

  py::class_<LmdbSink, IDynamicSink, std::shared_ptr<LmdbSink>>(md, "LmdbSink")
    .def(py::init<
           std::shared_ptr<lahuta::LMDBDatabase>, std::size_t
         >(),
         py::arg("db"), py::arg("batch_size") = 1024u);
}

} // namespace lahuta::bindings

#endif // LAHUTA_BINDINGS_DYNAMIC_SINKS_HPP
