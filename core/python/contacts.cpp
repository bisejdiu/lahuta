#include <optional>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "contacts.hpp"
#include "contacts/arpeggio/provider.hpp"
#include "contacts/engine.hpp"
#include "contacts/molstar/provider.hpp"

#include "entities/contact.hpp"
#include "entities/entity_id.hpp"
#include "entities/interaction_types.hpp"
#include "topology.hpp"

using namespace lahuta;

// clang-format off
void bind_contacts(py::module &m) {
  py::enum_<Category>(m, "Category")
    .value("None",                   Category::None)
    .value("None_",                  Category::None)
    .value("Generic",                Category::Generic)
    .value("Hydrophobic",            Category::Hydrophobic)
    .value("Halogen",                Category::Halogen)
    .value("HydrogenBond",           Category::HydrogenBond)
    .value("WeakHydrogenBond",       Category::WeakHydrogenBond)
    .value("PolarHydrogenBond",      Category::PolarHydrogenBond)
    .value("WeakPolarHydrogenBond",  Category::WeakPolarHydrogenBond)
    .value("Aromatic",               Category::Aromatic)
    .value("Ionic",                  Category::Ionic)
    .value("MetalCoordination",      Category::MetalCoordination)
    .value("CationPi",               Category::CationPi)
    .value("PiStacking",             Category::PiStacking)
    .value("Carbonyl",               Category::Carbonyl)
    .value("VanDerWaals",            Category::VanDerWaals)
    .value("DonorPi",                Category::DonorPi)
    .value("SulphurPi",              Category::SulphurPi)
    .value("CarbonPi",               Category::CarbonPi);

  py::enum_<Flavor>(m, "Flavor")
    .value("Default",  Flavor::Default)
    .value("Parallel", Flavor::Parallel)
    .value("TShape",   Flavor::TShape);

  py::class_<InteractionType> InteractionType_(m, "InteractionType");
  InteractionType_
    .def(py::init<Category, Flavor>(), py::arg("category") = Category::None, py::arg("flavor") = Flavor::Default)
    .def_readonly_static("None",                  &InteractionType::None)
    .def_readonly_static("None_",                 &InteractionType::None)
    .def_readonly_static("Generic",               &InteractionType::Generic)
    .def_readonly_static("Hydrophobic",           &InteractionType::Hydrophobic)
    .def_readonly_static("Halogen",               &InteractionType::Halogen)
    .def_readonly_static("Ionic",                 &InteractionType::Ionic)
    .def_readonly_static("CationPi",              &InteractionType::CationPi)
    .def_readonly_static("HydrogenBond",          &InteractionType::HydrogenBond)
    .def_readonly_static("WeakHydrogenBond",      &InteractionType::WeakHydrogenBond)
    .def_readonly_static("PolarHydrogenBond",     &InteractionType::PolarHydrogenBond)
    .def_readonly_static("WeakPolarHydrogenBond", &InteractionType::WeakPolarHydrogenBond)
    .def_readonly_static("MetalCoordination",     &InteractionType::MetalCoordination)
    .def_readonly_static("Aromatic",              &InteractionType::Aromatic)
    .def_readonly_static("PiStacking",            &InteractionType::PiStacking)
    .def_readonly_static("PiStackingP",           &InteractionType::PiStackingP)
    .def_readonly_static("PiStackingT",           &InteractionType::PiStackingT)
    .def_readonly_static("Carbonyl",              &InteractionType::Carbonyl)
    .def_readonly_static("VanDerWaals",           &InteractionType::VanDerWaals)
    .def_readonly_static("DonorPi",               &InteractionType::DonorPi)
    .def_readonly_static("SulphurPi",             &InteractionType::SulphurPi)
    .def_readonly_static("CarbonPi",              &InteractionType::CarbonPi)
    .def_property_readonly("category", [](const InteractionType &self){ return self.category; })
    .def_property_readonly("flavor",   [](const InteractionType &self){ return self.flavor; })
    .def("__repr__", [](const InteractionType &self){ return "InteractionType(" + interaction_type_to_string(self) + ")"; })
    .def("__str__",  [](const InteractionType &self){ return interaction_type_to_string(self); })
    .def("__eq__",   &InteractionType::operator==, py::is_operator())
    .def("__ne__",   &InteractionType::operator!=, py::is_operator());

  py::enum_<Kind>(m, "Kind")
    .value("Atom",  Kind::Atom)
    .value("Ring",  Kind::Ring)
    .value("Group", Kind::Group);

  py::class_<EntityID> EntityID_(m, "EntityID");
  EntityID_
    .def(py::init<uint64_t>())
    .def_property_readonly("kind",  &EntityID::kind)
    .def_property_readonly("index", &EntityID::index)
    .def("__int__", [](const EntityID &self)  { return EntityID::make(self.kind(), self.index()).raw; })
    .def("__repr__", [](const EntityID &self) { return "EntityID(" + self.to_string() + ")"; })
    .def("__str__",  [](const EntityID &self) { return self.to_string(); })
    .def(py::self == py::self)
    .def(py::self != py::self)
    .def(py::self <  py::self);

  py::class_<Contact> Contact_(m, "Contact");
  Contact_
    .def(py::init<EntityID, EntityID, float, InteractionType>(), py::arg("lhs"), py::arg("rhs"), py::arg("distance"), py::arg("type"))
    .def_readwrite("lhs",      &Contact::lhs)
    .def_readwrite("rhs",      &Contact::rhs)
    .def_readwrite("distance", &Contact::distance)
    .def_readwrite("type",     &Contact::type)
    .def("__eq__", &Contact::operator==, py::is_operator())
    .def("__ne__", &Contact::operator!=, py::is_operator())
    .def("__lt__", &Contact::operator<,  py::is_operator())
    .def("__repr__", [](const Contact &self) {
      return py::str("Contact(lhs={}, rhs={}, distance={}, type={})").format(self.lhs.to_string(), self.rhs.to_string(), self.distance, interaction_type_to_string(self.type));
    });

  py::class_<ContactSet> ContactSet_(m, "ContactSet");
  ContactSet_
    .def(py::init<>())
    .def("data",         [](ContactSet &self) { return self.data(); })
    .def("insert",       py::overload_cast<const ContactSet &>(&ContactSet::insert), py::arg("contacts"))
    .def("insert",       py::overload_cast<const Contact    &>(&ContactSet::insert), py::arg("contact" ))
    .def("size",         &ContactSet::size)
    .def("empty",        &ContactSet::empty)
    .def("make_generic", &ContactSet::make_generic)
    .def("__len__",      &ContactSet::size)
    .def("__and__",  [](const ContactSet &a, const ContactSet &b) { return a &  b; }, py::is_operator())
    .def("__or__",   [](const ContactSet &a, const ContactSet &b) { return a |  b; }, py::is_operator())
    .def("__sub__",  [](const ContactSet &a, const ContactSet &b) { return a -  b; }, py::is_operator())
    .def("__xor__",  [](const ContactSet &a, const ContactSet &b) { return a ^  b; }, py::is_operator())
    .def("__iand__", [](      ContactSet &a, const ContactSet &b) { return a &= b; }, py::is_operator())
    .def("__ior__",  [](      ContactSet &a, const ContactSet &b) { return a |= b; }, py::is_operator())
    .def("__isub__", [](      ContactSet &a, const ContactSet &b) { return a -= b; }, py::is_operator())
    .def("__ixor__", [](      ContactSet &a, const ContactSet &b) { return a ^= b; }, py::is_operator())
    .def("__getitem__", [](const ContactSet &self, int index) -> Contact {
      if (index < 0 || static_cast<size_t>(index) >= self.size()) throw py::index_error("Index out of range");
      return self.data()[static_cast<size_t>(index)];
    })
    .def("__iter__", [](const ContactSet &self) { return py::make_iterator(self.begin(), self.end()); }, py::keep_alive<0, 1>());

  // Engines
  using MsEngine = InteractionEngine<MolStarContactProvider>;
  using AgEngine = InteractionEngine<ArpeggioContactProvider>;

  py::class_<MsEngine>(m, "MolStarContactsEngine")
    .def(py::init<>())
    .def("compute", py::overload_cast<const Topology&>(&MsEngine::compute, py::const_), py::arg("topology"))
    .def("compute", py::overload_cast<const Topology&, std::optional<InteractionType>>(&MsEngine::compute, py::const_), py::arg("topology"), py::arg("only"));

  py::class_<AgEngine>(m, "ArpeggioContactsEngine")
    .def(py::init<>())
    .def("compute", py::overload_cast<const Topology&>(&AgEngine::compute, py::const_), py::arg("topology"))
    .def("compute", py::overload_cast<const Topology&, std::optional<InteractionType>>(&AgEngine::compute, py::const_), py::arg("topology"), py::arg("only"));
}
