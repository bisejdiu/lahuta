#include "lahuta.hpp"
#include "nsgrid.hpp"
#include <GraphMol/RDKitBase.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
// using namespace Lahuta;

void test_lahuta(py::module &_lahuta) {
  py::class_<Lahuta::GemmiSource> GemmiSource(_lahuta, "GemmiSource");
  py::class_<Lahuta::Luni> Luni(_lahuta, "Luni");

  GemmiSource
      .def(py::init<>())
      .def("process", &Lahuta::GemmiSource::process)
      .def("get_conformer", &Lahuta::GemmiSource::get_conformer);
  // .def("get_molecule",
  //      (RDKit::RWMol & (GemmiSource::*)()) & GemmiSource::get_molecule)
  // .def("get_molecule", (const RDKit::RWMol &(GemmiSource::*)() const) &
  //                          GemmiSource::get_molecule);

  Luni.def(py::init<std::string>())
      .def("get_neighbors", &Lahuta::Luni::get_neighbors)
      .def("find_neighbors", &Lahuta::Luni::find_neighbors)
      .def("get_cutoff", &Lahuta::Luni::get_cutoff);

  py::class_<NSResults>(_lahuta, "NSResults")
      .def(py::init<>())
      .def("add_neighbors", &NSResults::add_neighbors)
      .def("reserve_space", &NSResults::reserve_space)
      .def("get_neighbors", &NSResults::get_neighbors)
      .def("size", &NSResults::size)
      .def("filter", &NSResults::filter);
}


PYBIND11_MODULE(_lahuta, m) {
  m.doc() = "Lahuta: A Python binding for the Lahuta library";

  py::class_<RDKit::RWMol> LahutaRWMol(m, "RWMol");
  py::class_<RDKit::Conformer> LahutaConformer(m, "Conformer");

  test_lahuta(m);
}
