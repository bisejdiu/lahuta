#include "lahuta.hpp"
#include "nsgrid.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

// // Code to write Python bindings
// std::string file_name = "1a1e.cif.gz";
// Structure st = read_structure_gz(file_name);
// Lahuta::GemmiSource source = Lahuta::GemmiSource();
// source.process(st);
// Lahuta::Luni luni(source);

//
// py::class_<Lahuta::ISource>(m, "ISource")
//    .def("getMolecule", &Lahuta::ISource::getMolecule)
//    .def("getConformer", &Lahuta::ISource::getConformer)
//    .def("process", &Lahuta::ISource::process);

void test_lahuta(py::module &LT) {
  // ISource is a virtual class, so we can't instantiate it
  py::class_<Lahuta::GemmiSource> GemmiSource(LT, "GemmiSource");
  py::class_<Lahuta::Luni> Luni(LT, "Luni");

  // py::class_<Lahuta::GemmiSource, Lahuta::ISource>(LT, "GemmiSource")
  //     .def(py::init<>())
  //     .def("process", &Lahuta::GemmiSource::process);

  GemmiSource
      // .def(py::init<>(), py::arg("source") = Lahuta::GemmiSource())
      .def(py::init<>())
      .def("process", &Lahuta::GemmiSource::process);
      // .def("getMolecule", &Lahuta::GemmiSource::getMolecule)
      // .def_readonly("getMolecule", &Lahuta::GemmiSource::getMolecule)
      // .def("getConformer", &Lahuta::GemmiSource::getConformer);

  Luni
      // .def(py::init<Lahuta::GemmiSource>(), py::return_value_policy::copy)
      .def(py::init<std::string>())
      .def("getNeighborResults", &Lahuta::Luni::get_neighbors)
      .def("findNeighbors", &Lahuta::Luni::find_neighbors)
      .def("getCutoff", &Lahuta::Luni::get_cutoff);
      // .def("getMolecule", &Lahuta::Luni::getMolecule);

  py::class_<NSResults>(LT, "NSResults")
      .def(py::init<>())
      .def("addNeighbors", &NSResults::add_neighbors)
      .def("reserveSpace", &NSResults::reserve_space)
      .def("getNeighbors", &NSResults::get_neighbors)
      .def("getNeighborPairsSize", &NSResults::size)
      .def("filterByDistance", &NSResults::filter);
}


PYBIND11_MODULE(_lahuta, m) {
  m.doc() = "Lahuta: A Python binding for the Lahuta library";
  test_lahuta(m);
}

