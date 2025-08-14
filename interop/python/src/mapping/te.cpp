#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "mapper.hpp"

namespace py = pybind11;

// clang-format off
namespace lahuta::bindings {

void bind_te(py::module_ &m) {

  py::class_<TopologyMapper> topo_mapper(m, "TopologyMapper");

  py::enum_<TopologyMapper::MappingType>(topo_mapper, "MappingType")
    .value("Query",   TopologyMapper::MappingType::Query)
    .value("Target",  TopologyMapper::MappingType::Target);

  topo_mapper
  .def(py::init<SeqData&, TopologyMapper::MappingType>(), py::arg("sd"), py::arg("type"))
  .def("map", &TopologyMapper::map, py::arg("res"))
  .def("get_mapped_resid", &TopologyMapper::get_mapped_resid, py::arg("atom_index"))
  .def("get_entity_atoms", &TopologyMapper::get_entity_atoms, py::arg("entity"))
  .def("get_luni",         &TopologyMapper::get_luni, py::return_value_policy::reference_internal);

  py::class_<ContactEquivKey>(m, "ContactEquivKey")
    .def(py::init<>())
    .def_readwrite("e1_mapped_ids", &ContactEquivKey::e1_mapped_ids)
    .def_readwrite("e2_mapped_ids", &ContactEquivKey::e2_mapped_ids)
    .def_readwrite("contact_type",  &ContactEquivKey::contact_type);

  py::class_<EquivalencyConfig>(m, "EquivalencyConfig")
    .def(py::init<>())
    .def_readwrite("contact_resolution",  &EquivalencyConfig::contact_resolution)
    .def_readwrite("contact_type",        &EquivalencyConfig::contact_type)

    .def_readwrite("hbond_type",      &EquivalencyConfig::hbond_type)
    .def_readwrite("number_of_atoms", &EquivalencyConfig::number_of_atoms)
    .def_readwrite("atom_name",       &EquivalencyConfig::atom_name)
    .def_readwrite("element",         &EquivalencyConfig::element)
    .def_readwrite("resname",         &EquivalencyConfig::resname);

  py::class_<TopologicalEquivalency>(m, "TopologicalEquivalency")
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

  py::class_<LahutaMapper>(m, "LahutaMapper")
    .def(py::init<SeqData&, SeqData&>(),      py::arg("query"), py::arg("target"))
    .def("map", &LahutaMapper::map,           py::arg("res"),   py::arg("config") = std::nullopt)
    .def("evaluate", &LahutaMapper::evaluate, py::arg("c1"),    py::arg("m1"), py::arg("c2"), py::arg("m2"))
    ;
}
} // namespace lahuta::bindings
