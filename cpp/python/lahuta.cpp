#include "lahuta.hpp"
// #include "atom_types.hpp"
#include "nsgrid.hpp"
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <rdkit/GraphMol/RDKitBase.h>

namespace py = pybind11;
using namespace lahuta;

void test_lahuta(py::module &_lahuta) {
  py::class_<GemmiSource> GemmiSource(_lahuta, "GemmiSource");
  py::class_<Luni> Luni(_lahuta, "Luni");

  py::enum_<AtomType> AtomType(_lahuta, "AtomType");
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
      .def(pybind11::self | pybind11::self)
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
      .def_property_readonly(
          "value",
          [](enum AtomType flag) { return static_cast<uint32_t>(flag); })


      .def(~pybind11::self);

  AtomType.def("has", &has, "Check if flags contain a specific flag");
  AtomType.def("all", &all, "Check if flags contain all flags");
  AtomType.def("any", &any, "Check if flags contain any flags");
  AtomType.def("none", &none, "Check if flags contain none of the flags");
  AtomType.def("print_flags", &print_flags, "Print flags");
  AtomType
      .def("empty", &empty, "Check if flags are empty")
      // add __str__ method
      .def("__str__",
           [](enum AtomType flag) { return atom_type_to_string(flag); })
      .def("__repr__", [](enum AtomType flag) {
        return "<AtomType." + atom_type_to_string(flag) + ">";
      });

  // pybind11::implicitly_convertible<AtomType, uint32_t>();
  //
  // pybind11::class_<AtomType>(_lahuta, "AtomTypeOps")

  GemmiSource.def(py::init<>())
      .def("process", &GemmiSource::process)
      .def("get_conformer", &GemmiSource::get_conformer);
  // .def("get_molecule",
  //      (RDKit::RWMol & (GemmiSource::*)()) & GemmiSource::get_molecule)
  // .def("get_molecule", (const RDKit::RWMol &(GemmiSource::*)() const) &
  //                          GemmiSource::get_molecule);
  py::class_<NSResults>(_lahuta, "NSResults")
      .def(py::init<>())
      .def("add_neighbors", &NSResults::add_neighbors)
      .def("reserve_space", &NSResults::reserve_space)
      .def("size", &NSResults::size)
      .def("filter", &NSResults::filter)
      .def("type_filter", &NSResults::type_filter)
      .def("clear", &NSResults::clear)
      .def("get_pairs", &NSResults::get_pairs)
      .def("get_distances", &NSResults::get_distances)
      .def("get_luni", &NSResults::get_luni);
      // .def_property_readonly("luni", &NSResults::luni);

  Luni.def(py::init<std::string>())
      .def("find_neighbors", &Luni::find_neighbors)
      .def("get_atom_types", &Luni::get_atom_types)
      .def("match_smarts_string", &Luni::match_smarts_string)
      .def("get_cutoff", &Luni::get_cutoff);



    _lahuta.def("atom_type_to_string", &atom_type_to_string);
}

PYBIND11_MODULE(_lahuta, m) {
  m.doc() = "lahuta: A Python binding for the Lahuta library";

  py::class_<RDKit::RWMol> lahutaRWMol(m, "RWMol");
  py::class_<RDKit::Conformer> lahutaConformer(m, "Conformer");

  test_lahuta(m);
}
