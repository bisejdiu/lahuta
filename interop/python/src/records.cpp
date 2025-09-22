#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "entities/records.hpp"
#include "typing/types.hpp"

namespace py = pybind11;

// clang-format off
namespace lahuta::bindings {

void bind_records(py::module &m) {
  py::enum_<FeatureGroup>(m, "FeatureGroup", "Functional group classification used by the topology")
    .value("None_",           FeatureGroup::None)
    .value("QuaternaryAmine", FeatureGroup::QuaternaryAmine)
    .value("TertiaryAmine",   FeatureGroup::TertiaryAmine)
    .value("Sulfonium",       FeatureGroup::Sulfonium)
    .value("SulfonicAcid",    FeatureGroup::SulfonicAcid)
    .value("Sulfate",         FeatureGroup::Sulfate)
    .value("Phosphate",       FeatureGroup::Phosphate)
    .value("Halocarbon",      FeatureGroup::Halocarbon)
    .value("Guanidine",       FeatureGroup::Guanidine)
    .value("Acetamidine",     FeatureGroup::Acetamidine)
    .value("Carboxylate",     FeatureGroup::Carboxylate);

  py::class_<AtomRec>(m, "AtomRec", "Atom record with typing and RDKit atom handle")
    .def_property("type",
      [](const AtomRec &self) { return self.type; },
      [](AtomRec &self, const AtomType &value) { self.type = value; },
      "Atom type classification"
    )
    .def("idx", [](const AtomRec &self) { return self.atom.get().getIdx(); }, "RDKit atom index")
    .def("__repr__", [](const AtomRec &self) {
        return py::str("AtomRec(type={}, idx={})").format(static_cast<uint32_t>(self.type), self.atom.get().getIdx() );
    });

  py::class_<RingRec>(m, "RingRec", "Ring record with atom membership")
    .def(py::init<>(), "Create empty ring record")
    .def_readwrite("atoms",    &RingRec::atoms,    "Atoms participating in the ring")
    .def_readwrite("aromatic", &RingRec::aromatic, "Whether the ring is aromatic")
    .def_property_readonly("size", [](const RingRec &self) { return static_cast<std::size_t>(self.atoms.size()); }, "Number of atoms in the ring")
    .def("__repr__", [](const RingRec &self) {
        return py::str("RingRec(atoms={}, aromatic={})").format(self.atoms.size(), self.aromatic);
    });

  py::class_<GroupRec>(m, "GroupRec", "Functional group record with atoms")
    .def(py::init<>(), "Create empty group record")
    .def_property(
      "a_type", // FIX: a_type vs type naming inconsistency
      [](const GroupRec &self) { return self.a_type; },
      [](GroupRec &self, const AtomType &value) { self.a_type = value; },
      "Atom type associated with this group"
    )
    .def_property(
      "type",
      [](const GroupRec &self) { return self.type; },
      [](GroupRec &self, const FeatureGroup &value) { self.type = value; },
      "Functional group classification"
    )
    .def_readwrite("atoms",  &GroupRec::atoms,  "Atoms participating in the group")
    .def("__repr__", [](const GroupRec &self) {
        return py::str("GroupRec(a_type={}, type={}, atoms={})").format(
            static_cast<uint32_t>(self.a_type),
            static_cast<int>(self.type),
            self.atoms.size()
        );
    });
}

} // namespace lahuta::bindings
