#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "fseek/ops.hpp"
#include "fseek/prefilter.hpp"
#include "numpy_utils.hpp"
#include "procs/search_and_align.hpp"

// clang-format off
namespace py = pybind11;

namespace lahuta::bindings {

void bind_align(py::module_ &m) {
    py::class_<LahutaAlignerBase,  std::shared_ptr<LahutaAlignerBase>>  a_base(m, "LahutaAlignerBase");

    py::class_<ProcessingConfig> (m, "ProcessingConfig")
        .def(py::init<>())
        .def_readwrite("query_chunk_size",  &ProcessingConfig::query_chunk_size)
        .def_readwrite("target_chunk_size", &ProcessingConfig::target_chunk_size)
        .def_readwrite("allow_self_ops",    &ProcessingConfig::allow_self_ops);

    py::class_<LahutaProcessor,    std::shared_ptr<LahutaProcessor>>(m, "LahutaProcessor")
        .def(py::init<std::shared_ptr<LahutaAlignerBase>, const ProcessingConfig &>(), py::arg("aligner"), py::arg("config"))
        .def("process",
             py::overload_cast<const std::vector<std::string>&, const std::vector<std::string>&>(&LahutaProcessor::process),
             py::arg("query_files"), py::arg("target_files"))
        .def("process",
             py::overload_cast<const std::vector<std::string>&>(&LahutaProcessor::process),
             py::arg("query_files"));

    py::class_<SeqData, std::shared_ptr<SeqData>>(m, "SeqData")
        .def(py::init<>())
        .def("size", (int (SeqData::*)() const) &SeqData::size, "Returns the number of residues in the sequence")
        .def_readwrite("Seq3Di",      &SeqData::Seq3Di)
        .def_readwrite("SeqAA",       &SeqData::SeqAA)
        .def_readwrite("CaData",      &SeqData::CaData)
        .def_readwrite("file_name",   &SeqData::file_name)
        .def_readwrite("chain_name",  &SeqData::chain_name)
        .def_property_readonly("x", [](SeqData &self) {
            size_t n = self.SeqAA.size();
            auto vec = std::vector<float>(self.CaData.begin(), self.CaData.begin() + n);
            return numpy::as_numpy_copy(vec);
        })
        .def_property_readonly("y", [](SeqData &self) {
            size_t n = self.SeqAA.size();
            auto vec = std::vector<float>(self.CaData.begin() + n, self.CaData.begin() + 2 * n);
            return numpy::as_numpy_copy(vec);
        })
        .def_property_readonly("z", [](SeqData &self) {
            size_t n = self.SeqAA.size();
            auto vec = std::vector<float>(self.CaData.begin() + 2 * n, self.CaData.begin() + 3 * n);
            return numpy::as_numpy_copy(vec);
        })
        .def("__lt__", [](const SeqData &a, const SeqData &b) { return a < b; });

    py::class_<AlignerResults, std::shared_ptr<AlignerResults>>(m, "AlignerResults")
        .def(py::init<>())
        .def_property("query",
             [](AlignerResults &self) { return self.query; },
             [](AlignerResults &self, std::shared_ptr<SeqData> q) { self.query = q; },
             py::keep_alive<1, 2>()) // Keep AlignerResults alive as long as returned query is alive
        .def_property("target",
             [](AlignerResults &self) { return self.target; },
             [](AlignerResults &self, std::shared_ptr<SeqData> t) { self.target = t; },
             py::keep_alive<1, 2>()) // Keep AlignerResults alive as long as returned target is alive
        .def_readwrite("results", &AlignerResults::results);

    py::class_<LahutaAligner, LahutaAlignerBase,  std::shared_ptr<LahutaAligner>>(m, "LahutaAligner")
        .def(py::init<FoldSeekOps, PrefilterOptions, unsigned int>(),
             py::arg_v("ops",    FoldSeekOps(),      "FoldSeekOps()"),
             py::arg_v("pf_ops", PrefilterOptions(), "PrefilterOptions()"),
             py::arg("n_threads") = 0)
        .def("run",         &LahutaAligner::run, py::arg("query_files"), py::arg("target_files"))
        .def("get_results", &LahutaAligner::get_results);
}

} // namespace lahuta::bindings
