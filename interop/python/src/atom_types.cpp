/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   constexpr char p1[] = "besian", p2[] = "sejdiu", p3[] = "@gmail.com"; std::string s;
 *   s.append(std::begin(p1), std::end(p1) - 1);
 *   s.append(std::begin(p2), std::end(p2) - 1);
 *   s.append(std::begin(p3), std::end(p3) - 1);
 *   return s;
 * }();
 *
 */

#include <pybind11/attr.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "typing/flags.hpp"
#include "typing/types.hpp"

namespace py = pybind11;

// clang-format off
namespace lahuta::bindings {
namespace Flags = AtomTypeFlags;

void bind_atom_types(py::module &m) {
  py::enum_<AtomType> (m, "AtomType")
    .value("None_",             AtomType::None)
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
    .value("Invalid",           AtomType::Invalid)

    .def("__or__",  [](AtomType  lhs, AtomType rhs) { return lhs |  rhs; }, py::is_operator())
    .def("__and__", [](AtomType  lhs, AtomType rhs) { return lhs &  rhs; }, py::is_operator())
    .def("__xor__", [](AtomType  lhs, AtomType rhs) { return lhs ^  rhs; }, py::is_operator())
    .def("__ior__", [](AtomType &lhs, AtomType rhs) { return lhs |= rhs; }, py::is_operator())
    .def("__iand__",[](AtomType &lhs, AtomType rhs) { return lhs &= rhs; }, py::is_operator())
    .def("__ixor__",[](AtomType &lhs, AtomType rhs) { return lhs ^= rhs; }, py::is_operator())

    .def_property_readonly("label", &lahuta::atom_type_to_string, "Human-readable label for this atom type")
    .def("has",     [](AtomType self, AtomType flag)  { return Flags::has(self, flag); },      "Return True if all bits in flag are set in self",  py::arg("flag"))
    .def("has_any", [](AtomType self, AtomType flags) { return Flags::has_any(self, flags); }, "Return True if any bits in flags are set in self", py::arg("flags"))
    .def("all",     [](AtomType self, AtomType flags) { return Flags::all(self, flags); },     "Return True if self contains all bits in flags",   py::arg("flags"))
    .def("any",     [](AtomType self, AtomType flags) { return Flags::any(self, flags); },     "Return True if self contains any bits in flags",   py::arg("flags"))
    .def("none",    [](AtomType self, AtomType flags) { return Flags::none(self, flags); },    "Return True if self contains none of the bits in flags", py::arg("flags"))
    .def_property_readonly("empty",   [](AtomType self) { return Flags::empty(self); },        "True if self has no bits set (AtomType.None_)")
    .def("split",   [](AtomType self) { return Flags::split(self); },                          "Return list of individual flags set in self")
    ;
}
} // namespace lahuta::bindings
