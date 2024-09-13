#include "atom_types.hpp"
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <rdkit/Geometry/point.h>
#include <rdkit/GraphMol/RDKitBase.h>

namespace py = pybind11;
using namespace lahuta;
using namespace AtomTypeFlags;

void xbind_atom_types(py::module &m) {

  py::enum_<AtomType> AtomType(m, "AtomType");
  AtomType.value("NONE", AtomType::NONE)
      .value("HBOND_ACCEPTOR", AtomType::HBOND_ACCEPTOR)
      .value("HBOND_DONOR", AtomType::HBOND_DONOR)
      .value("WEAK_HBOND_ACCEPTOR", AtomType::WEAK_HBOND_ACCEPTOR)
      .value("WEAK_HBOND_DONOR", AtomType::WEAK_HBOND_DONOR)
      .value("POS_IONISABLE", AtomType::POS_IONISABLE)
      .value("NEG_IONISABLE", AtomType::NEG_IONISABLE)
      .value("CARBONYL_OXYGEN", AtomType::CARBONYL_OXYGEN)
      .value("CARBONYL_CARBON", AtomType::CARBONYL_CARBON)
      .value("AROMATIC", AtomType::AROMATIC)
      .value("HYDROPHOBIC", AtomType::HYDROPHOBIC)
      .value("XBOND_ACCEPTOR", AtomType::XBOND_ACCEPTOR)
      .value("XBOND_DONOR", AtomType::XBOND_DONOR)
      .value("INVALID", AtomType::INVALID)
      .export_values();

  AtomType.def(pybind11::self | pybind11::self)
      .def(pybind11::self & pybind11::self)
      .def(pybind11::self ^ pybind11::self)
      .def(pybind11::self |= pybind11::self)
      .def(pybind11::self &= pybind11::self)
      .def(pybind11::self ^= pybind11::self)
      .def_property_readonly(
          "name", [](enum AtomType flag) { return atom_type_to_string(flag); })
      .def("components",
           [](enum AtomType flag) {
             std::vector<enum AtomType> components;
             if (has(flag, AtomType::HBOND_ACCEPTOR))
               components.push_back(AtomType::HBOND_ACCEPTOR);
             if (has(flag, AtomType::HBOND_DONOR))
               components.push_back(AtomType::HBOND_DONOR);
             if (has(flag, AtomType::WEAK_HBOND_ACCEPTOR))
               components.push_back(AtomType::WEAK_HBOND_ACCEPTOR);
             if (has(flag, AtomType::WEAK_HBOND_DONOR))
               components.push_back(AtomType::WEAK_HBOND_DONOR);
             if (has(flag, AtomType::POS_IONISABLE))
               components.push_back(AtomType::POS_IONISABLE);
             if (has(flag, AtomType::NEG_IONISABLE))
               components.push_back(AtomType::NEG_IONISABLE);
             if (has(flag, AtomType::CARBONYL_OXYGEN))
               components.push_back(AtomType::CARBONYL_OXYGEN);
             if (has(flag, AtomType::CARBONYL_CARBON))
               components.push_back(AtomType::CARBONYL_CARBON);
             if (has(flag, AtomType::AROMATIC))
               components.push_back(AtomType::AROMATIC);
             if (has(flag, AtomType::HYDROPHOBIC))
               components.push_back(AtomType::HYDROPHOBIC);
             if (has(flag, AtomType::XBOND_ACCEPTOR))
               components.push_back(AtomType::XBOND_ACCEPTOR);
             if (has(flag, AtomType::XBOND_DONOR))
               components.push_back(AtomType::XBOND_DONOR);
             if (has(flag, AtomType::INVALID))
               components.push_back(AtomType::INVALID);
             return components;
           })
      .def_property_readonly("value", [](enum AtomType flag) { 
        auto flag_str = atom_type_to_string(flag);
        if (flag_str.back() == ' ') {
          flag_str.pop_back();
        }
        return flag_str;
        // return static_cast<uint32_t>(flag); 
      })

      .def(~pybind11::self);

  AtomType.def("has", &has, "Check if flags contain a specific flag");
  AtomType.def("has_enum_as_str", &has_enum_as_string,
               "Check if flags contain a specific flag");
  AtomType.def("get_enum_as_str", &get_enum_as_string, "Get a specific flag");
  AtomType.def("all", &all, "Check if flags contain all flags");
  AtomType.def("any", &any, "Check if flags contain any flags");
  AtomType.def("none", &none, "Check if flags contain none of the flags");
  AtomType.def("print_flags", &print_flags, "Print flags");
  AtomType.def("empty", &empty, "Check if flags are empty")
      .def("__str__",
           [](enum AtomType flag) { return atom_type_to_string(flag); })
      .def("__repr__", [](enum AtomType flag) {
        return "<AtomType." + atom_type_to_string(flag) + ">";
      });

  // m.def("get_enum_as_str", &get_enum_as_string, "Get a specific flag");
  m.def("atom_type_to_string", &atom_type_to_string);
}
