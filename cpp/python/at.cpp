#include <pybind11/attr.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "at.hpp"
#include "atom_types.hpp"

using namespace lahuta;
namespace Flags = AtomTypeFlags;
namespace py = pybind11;

// clang-format off

void bind_atom_types(py::module &m) {
  auto ATFlags = m.def_submodule("Flags", "Bindings for common utilities");
  py::enum_<AtomType> AtomType(m, "AtomType");

  AtomType
    .value("NONE",                AtomType::NONE)
    .value("HBOND_ACCEPTOR",      AtomType::HBOND_ACCEPTOR)
    .value("HBOND_DONOR",         AtomType::HBOND_DONOR)
    .value("WEAK_HBOND_ACCEPTOR", AtomType::WEAK_HBOND_ACCEPTOR)
    .value("WEAK_HBOND_DONOR",    AtomType::WEAK_HBOND_DONOR)
    .value("POS_IONISABLE",       AtomType::POS_IONISABLE)
    .value("NEG_IONISABLE",       AtomType::NEG_IONISABLE)
    .value("CARBONYL_OXYGEN",     AtomType::CARBONYL_OXYGEN)
    .value("CARBONYL_CARBON",     AtomType::CARBONYL_CARBON)
    .value("AROMATIC",            AtomType::AROMATIC)
    .value("HYDROPHOBIC",         AtomType::HYDROPHOBIC)
    .value("XBOND_ACCEPTOR",      AtomType::XBOND_ACCEPTOR)
    .value("XBOND_DONOR",         AtomType::XBOND_DONOR)
    .value("INVALID",             AtomType::INVALID);

  AtomType
    .def(py::self |  py::self)
    .def(py::self &  py::self)
    .def(py::self ^  py::self)
    .def(py::self |= py::self)
    .def(py::self &= py::self)
    .def(py::self ^= py::self)
    .def(~py::self)
    .def_property_readonly("value", [](const enum AtomType flag) { return atom_type_to_string(flag); })
    .def("split",                   [](const enum AtomType flag) { return Flags::split(flag); })
    // See: https://github.com/pybind/pybind11/issues/2537#issuecomment-702941967
    .def("__str__", [](enum AtomType flag) { return "AtomType(" + atom_type_to_string(flag) + ")"; }, py::prepend());

  ATFlags
    .def("has",     &Flags::has,     "Check if flags contain a specific flag",   py::arg("types"), py::arg("type"))
    .def("has_any", &Flags::has_any, "Check if flags contain any of the flags",  py::arg("types"), py::arg("type"))
    .def("all",     &Flags::all,     "Check if flags contain all flags",         py::arg("types"), py::arg("type"))
    .def("any",     &Flags::any,     "Check if flags contain any flags",         py::arg("types"), py::arg("type"))
    .def("none",    &Flags::none,    "Check if flags contain none of the flags", py::arg("types"), py::arg("type"))
    .def("empty",   &Flags::empty,   "Check if flags are empty",                 py::arg("types"))
    .def("split",   &Flags::split,   "Split flags into components",              py::arg("types"))
    .def("string_to_atom_type", &string_to_atom_type, "Convert string to atom type", py::arg("types"))
    .def("get_enum_as_string",  &get_enum_as_string,  "Get enum as string",          py::arg("types"));

}
