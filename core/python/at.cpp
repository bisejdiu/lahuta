#include <pybind11/attr.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "at.hpp"
#include "typing/types.hpp"
#include "typing/flags.hpp"

using namespace lahuta;
namespace Flags = AtomTypeFlags;
namespace py = pybind11;

// clang-format off

inline AtomType string_to_atom_type(const std::string &flag_name) {
  static const std::unordered_map<std::string, AtomType> stringToEnum = {
      {"None",              AtomType::None},
      {"HbondAcceptor",     AtomType::HbondAcceptor},
      {"HbondDonor",        AtomType::HbondDonor},
      {"WeakHbondAcceptor", AtomType::WeakHbondAcceptor},
      {"WeakHbondDonor",    AtomType::WeakHbondDonor},
      {"PositiveCharge",    AtomType::PositiveCharge},
      {"NegativeCharge",    AtomType::NegativeCharge},
      {"CarbonylOxygen",    AtomType::CarbonylOxygen},
      {"CarbonylCarbon",    AtomType::CarbonylCarbon},
      {"Aromatic",          AtomType::Aromatic},
      {"Hydrophobic",       AtomType::Hydrophobic},
      {"XBondAcceptor",     AtomType::XBondAcceptor},
      {"XbondDonor",        AtomType::XbondDonor},
      {"IonicTypePartner",  AtomType::IonicTypePartner},
      {"DativeBondPartner", AtomType::DativeBondPartner},
      {"TransitionMetal",   AtomType::TransitionMetal},
      {"IonicTypeMetal",    AtomType::IonicTypeMetal},
      {"Invalid",           AtomType::Invalid}};

  auto it = stringToEnum.find(flag_name);
  if (it != stringToEnum.end()) {
    return it->second;
  }
  throw std::invalid_argument("Invalid AtomType flag name: " + flag_name);
}

inline AtomType get_enum_as_string(std::string flag_name) { return string_to_atom_type(flag_name); }


void bind_atom_types(py::module &m) {
  auto ATFlags = m.def_submodule("Flags", "Bindings for common utilities");
  py::enum_<AtomType> AtomType(m, "AtomType");

  AtomType
    .value("None",              AtomType::None)
    .value("HbondAcceptor",     AtomType::HbondAcceptor)
    .value("HbondDonor",        AtomType::HbondDonor)
    .value("WeakHbondAcceptor", AtomType::WeakHbondAcceptor)
    .value("WeakHbondDonor",    AtomType::WeakHbondDonor)
    .value("PositiveCharge",    AtomType::PositiveCharge)
    .value("NegativeCharge",    AtomType::NegativeCharge)
    .value("CarbonylOxygen",    AtomType::CarbonylOxygen)
    .value("CarbonylCarbon",    AtomType::CarbonylCarbon)
    .value("Aromatic",          AtomType::Aromatic)
    .value("Hydrophobic",       AtomType::Hydrophobic)
    .value("XBondAcceptor",     AtomType::XBondAcceptor)
    .value("XbondDonor",        AtomType::XbondDonor)
    .value("IonicTypePartner",  AtomType::IonicTypePartner)
    .value("DativeBondPartner", AtomType::DativeBondPartner)
    .value("TransitionMetal",   AtomType::TransitionMetal)
    .value("IonicTypeMetal",    AtomType::IonicTypeMetal)
    .value("Invalid",           AtomType::Invalid);


  AtomType
    .def(py::self |  py::self)
    .def(py::self &  py::self)
    .def(py::self ^  py::self)
    .def(py::self |= py::self)
    .def(py::self &= py::self)
    .def(py::self ^= py::self)
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
