#include "contacts.hpp"
#include "contacts/interactions.hpp"
#include "entities/entity_id.hpp"
#include "entities/contact.hpp"
#include "topology.hpp"

using namespace lahuta;

// clang-format off
void bind_contacts(py::module &m) {
  py::class_<InteractionOptions> InteractionOptions_(m, "InteractionOptions");
  py::class_<Interactions>       Interactions_      (m, "Interactions");
  py::class_<InteractionType>     InteractionType_   (m, "InteractionType");

  py::enum_<Kind> EntityType_(m, "EntityType");
  py::class_<EntityID>    EntityID_  (m, "EntityID");
  py::class_<Contact>     Contact_   (m, "Contact");
  py::class_<ContactSet>  Contacts_  (m, "Contacts");


  InteractionOptions_
    .def(py::init<float>(), py::arg("distance_cutoff") = 5.0f)

    .def_readwrite("hbond",       &InteractionOptions::hbond)
    .def_readwrite("weak_hbond",  &InteractionOptions::weak_hbond)
    .def_readwrite("hydrophobic", &InteractionOptions::hydrophobic)
    .def_readwrite("ionic",       &InteractionOptions::ionic)
    .def_readwrite("metalic",     &InteractionOptions::metalic)
    .def_readwrite("cationpi",    &InteractionOptions::cationpi)
    .def_readwrite("pistacking",  &InteractionOptions::pistacking)
    .def_readwrite("halogen",     &InteractionOptions::halogen);
    // .def_readwrite("sort_globally",         &InteractionOptions::sort_globally);


  Interactions_
    .def(py::init<const class Topology &, InteractionOptions>(), py::arg("topology"), py::arg("opts") = InteractionOptions{})

    .def("hbond",       &Interactions::hbond)
    .def("weak_hbond",  &Interactions::weak_hbond)
    .def("hydrophobic", &Interactions::hydrophobic)
    .def("halogen",     &Interactions::halogen)
    .def("ionic",       &Interactions::ionic)
    .def("metalic",     &Interactions::metalic)
    .def("cationpi",    &Interactions::cationpi)
    .def("pistacking",  &Interactions::pistacking)

    .def("compute_contacts", &Interactions::compute_contacts);

    InteractionType_
      .def_readonly_static("None",             &InteractionType::None)
      .def_readonly_static("Hydrophobic",      &InteractionType::Hydrophobic)
      .def_readonly_static("HydrogenBond",     &InteractionType::HydrogenBond)
      .def_readonly_static("WeakHydrogenBond", &InteractionType::WeakHydrogenBond)
      .def_readonly_static("Ionic",            &InteractionType::Ionic)
      .def_readonly_static("MetalCoordination",&InteractionType::MetalCoordination)
      .def_readonly_static("CationPi",         &InteractionType::CationPi)
      .def_readonly_static("PiStackingP",      &InteractionType::PiStackingP)
      .def_readonly_static("PiStackingT",      &InteractionType::PiStackingT);


    EntityType_
      .value("Atom",  Kind::Atom)
      .value("Ring",  Kind::Ring)
      .value("Group", Kind::Group);


    EntityID_ 
      .def(py::init<uint64_t>())
      .def("__int__", [](const EntityID &self) { return EntityID::make(self.kind(), self.index()); });


    Contact_
      .def(py::init<EntityID, EntityID, float, InteractionType>(), py::arg("index1"), py::arg("index2"), py::arg("distance"), py::arg("type"))

      .def_readwrite("index1",   &Contact::lhs)
      .def_readwrite("index2",   &Contact::rhs)
      .def_readwrite("distance", &Contact::distance)
      .def_readwrite("type",     &Contact::type)

      .def("__eq__", &Contact::operator==, py::is_operator())
      .def("__ne__", &Contact::operator!=, py::is_operator())
      .def("__lt__", &Contact::operator<,  py::is_operator())
      // .def("__gt__", &Contact::operator>,  py::is_operator())
      // .def("__le__", &Contact::operator<=, py::is_operator())
      // .def("__ge__", &Contact::operator>=, py::is_operator())

      .def("__str__", [](const Contact &self) {
        return "<Contact index1=" + self.lhs.to_string() //
               + ", index2=" + self.rhs.to_string() //
               + ", distance=" + std::to_string(self.distance)
               + ", type=" + std::to_string(static_cast<int>(self.type)) + ">";
      });


    Contacts_
      .def(py::init<>())
      // .def(py::init<const class Luni *>(), py::arg("luni"))

      .def("interactions",  [](ContactSet &self) { return self.data(); })
      .def("is_sorted",     [](ContactSet &self) { return true; })

      .def("add", py::overload_cast<const ContactSet &>(&ContactSet::insert), py::arg("contacts"))
      .def("add", py::overload_cast<const Contact &> (&ContactSet::insert), py::arg("interaction"))

      .def("size",         &ContactSet::size)
      .def("make_generic", &ContactSet::make_generic)

      // .def("__eq__",   &ContactSet::operator==, py::is_operator())
      // .def("__ne__",   &ContactSet::operator!=, py::is_operator())
      .def("__and__",  [](ContactSet &self, ContactSet &other) {return self &  other;}, py::is_operator())
      .def("__or__",   [](ContactSet &self, ContactSet &other) {return self |  other;}, py::is_operator())
      .def("__sub__",  [](ContactSet &self, ContactSet &other) {return self -  other;}, py::is_operator())
      .def("__xor__",  [](ContactSet &self, ContactSet &other) {return self ^  other;}, py::is_operator())
      .def("__iand__", [](ContactSet &self, ContactSet &other) {return self &= other;}, py::is_operator())
      .def("__ior__",  [](ContactSet &self, ContactSet &other) {return self |= other;}, py::is_operator())
      .def("__isub__", [](ContactSet &self, ContactSet &other) {return self -= other;}, py::is_operator())
      .def("__ixor__", [](ContactSet &self, ContactSet &other) {return self ^= other;}, py::is_operator())

       // FIX: add support for slicing
      .def("__getitem__", [](const ContactSet &self, int index) -> Contact {
        if (index < 0 || index >= self.size()) { throw py::index_error("Index out of range"); }
            return self.data()[index];
       });
}
