#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "array.hpp"
#include "mapper.hpp"
#include "py_te.hpp"

namespace py = pybind11;
using namespace lahuta;

// clang-format off
void bind_te(py::module &_lahuta) {

  py::class_<LahutaMapper>    lahuta_mapper   (_lahuta, "LahutaMapper");
  py::class_<TopologyMapper>  topology_mapper (_lahuta, "TopologyMapper");

  py::class_<ContactEquivKey>         contact_equivkey  (_lahuta, "ContactEquivKey");
  py::class_<EquivalencyConfig>       equivalency_config(_lahuta, "EquivalencyConfig");
  py::class_<TopologicalEquivalency>  topological_equiv (_lahuta, "TopologicalEquivalency");

  py::enum_<TopologyMapper::MappingType> mapping_type(topology_mapper, "MappingType");

  mapping_type
    .value("Query",   TopologyMapper::MappingType::Query)
    .value("Target",  TopologyMapper::MappingType::Target);

  topology_mapper
  .def(py::init<SeqData&, TopologyMapper::MappingType>(), py::arg("sd"), py::arg("type"))

  .def("map", &TopologyMapper::map, py::arg("res"))
  .def("get_mapped_resid", &TopologyMapper::get_mapped_resid, py::arg("atom_index"))
  .def("get_entity_atoms", &TopologyMapper::get_entity_atoms, py::arg("entity"))
  .def("get_luni",         &TopologyMapper::get_luni, py::return_value_policy::reference_internal);

  contact_equivkey
    .def(py::init<>())
    .def_readwrite("e1_mapped_ids", &ContactEquivKey::e1_mapped_ids)
    .def_readwrite("e2_mapped_ids", &ContactEquivKey::e2_mapped_ids)
    .def_readwrite("contact_type",  &ContactEquivKey::contact_type);

  equivalency_config
    .def(py::init<>())
    .def_readwrite("contact_resolution",  &EquivalencyConfig::contact_resolution)
    .def_readwrite("contact_type",        &EquivalencyConfig::contact_type)

    .def_readwrite("hbond_type",      &EquivalencyConfig::hbond_type)
    .def_readwrite("number_of_atoms", &EquivalencyConfig::number_of_atoms)
    .def_readwrite("atom_name",       &EquivalencyConfig::atom_name)
    .def_readwrite("element",         &EquivalencyConfig::element)
    .def_readwrite("resname",         &EquivalencyConfig::resname);

  topological_equiv
    .def(py::init<const TopologyMapper&, const TopologyMapper&, std::optional<EquivalencyConfig>>(),
         py::arg("lm1"), py::arg("lm2"), py::arg("config") = std::nullopt)
    .def("evaluate",
         py::overload_cast<const RDKit::Atom*, const RDKit::Atom*>(&TopologicalEquivalency::evaluate, py::const_),
         py::arg("a1"), py::arg("a2"))
    .def("evaluate_atoms",
         py::overload_cast<const std::vector<const RDKit::Atom *>, const std::vector<const RDKit::Atom *>>(&TopologicalEquivalency::evaluate, py::const_),
         py::arg("a1"), py::arg("a2"))
    .def("should_consider", &TopologicalEquivalency::should_consider, py::arg("atoms"), py::arg("mt"))
    .def("is_mappable",     &TopologicalEquivalency::is_mappable,     py::arg("mt"), py::arg("atom"))
    .def("get_lm",          &TopologicalEquivalency::get_lm,          py::arg("type"), py::return_value_policy::reference_internal);

  lahuta_mapper
    .def(py::init<SeqData&, SeqData&>(),      py::arg("query"), py::arg("target"))
    .def("map", &LahutaMapper::map,           py::arg("res"),   py::arg("config") = std::nullopt)
    .def("evaluate", &LahutaMapper::evaluate, py::arg("c1"),    py::arg("m1"), py::arg("c2"), py::arg("m2"))
    ;
}
