#include "contacts.hpp"
#include "contacts/interactions.hpp"
#include "lahuta.hpp"
#include "neighbors.hpp"

using namespace lahuta;
// clang-format off

void bind_contacts(py::module &m) {
  py::class_<InteractionOptions> InteractionOptions_(m, "InteractionOptions");
  py::class_<Interactions>       Interactions_      (m, "Interactions");
  py::enum_<InteractionType>     InteractionType_   (m, "InteractionType");

  py::enum_<EntityType> EntityType_(m, "EntityType");
  py::class_<EntityID>  EntityID_  (m, "EntityID");
  py::class_<Contact>   Contact_   (m, "Contact");
  py::class_<Contacts>  Contacts_  (m, "Contacts");


  InteractionOptions_
    .def(py::init<float>(), py::arg("distance_cutoff") = 5.0f)

    .def_readwrite("hbond",       &InteractionOptions::hbond)
    .def_readwrite("weak_hbond",  &InteractionOptions::weak_hbond)
    .def_readwrite("hydrophobic", &InteractionOptions::hydrophobic)
    .def_readwrite("ionic",       &InteractionOptions::ionic)
    .def_readwrite("metalic",     &InteractionOptions::metalic)
    .def_readwrite("cationpi",    &InteractionOptions::cationpi)
    .def_readwrite("pistacking",  &InteractionOptions::pistacking)
    .def_readwrite("halogen",     &InteractionOptions::halogen)

    .def_readwrite("sort_globally",         &InteractionOptions::sort_globally)
    .def_readwrite("sort_per_contact_type", &InteractionOptions::sort_per_contact_type);


  Interactions_
    .def(py::init<const class Luni &, std::optional<InteractionOptions>>(), py::arg("luni"), py::arg("opts") = std::nullopt)

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
      .value("None",             InteractionType::None)
      .value("Any",              InteractionType::Any)
      .value("Hydrophobic",      InteractionType::Hydrophobic)
      .value("HydrogenBond",     InteractionType::HydrogenBond)
      .value("WeakHydrogenBond", InteractionType::WeakHydrogenBond)
      .value("Ionic",            InteractionType::Ionic)
      .value("MetalCoordination",InteractionType::MetalCoordination)
      .value("CationPi",         InteractionType::CationPi)
      .value("PiStackingP",      InteractionType::PiStackingP)
      .value("PiStackingT",      InteractionType::PiStackingT);


    EntityType_
      .value("Atom",  EntityType::Atom)
      .value("Ring",  EntityType::Ring)
      .value("Group", EntityType::Group);


    EntityID_ 
      .def(py::init<uint64_t>())
      .def("__int__", [](const EntityID &self) { return static_cast<uint64_t>(self); });


    Contact_
      .def(py::init<EntityID, EntityID, float>(),                  py::arg("index1"), py::arg("index2"), py::arg("distance"))
      .def(py::init<EntityID, EntityID, float, InteractionType>(), py::arg("index1"), py::arg("index2"), py::arg("distance"), py::arg("type"))

      .def_readwrite("index1",   &Contact::entity1)
      .def_readwrite("index2",   &Contact::entity2)
      .def_readwrite("distance", &Contact::distance)
      .def_readwrite("type",     &Contact::type)

      .def("__eq__", &Contact::operator==, py::is_operator())
      .def("__ne__", &Contact::operator!=, py::is_operator())
      .def("__lt__", &Contact::operator<,  py::is_operator())
      .def("__gt__", &Contact::operator>,  py::is_operator())
      .def("__le__", &Contact::operator<=, py::is_operator())
      .def("__ge__", &Contact::operator>=, py::is_operator())

      .def("__str__", [](const Contact &self) {
        return "<Contact index1=" + std::to_string(self.entity1) //
               + ", index2=" + std::to_string(self.entity2)      //
               + ", distance=" + std::to_string(self.distance)
               + ", type=" + std::to_string(static_cast<int>(self.type)) + ">";
      });


    Contacts_
      .def(py::init<>())
      .def(py::init<const class Luni *>(), py::arg("luni"))

      .def_readwrite("interactions",  &Contacts::interactions)
      .def_readwrite("is_sorted",     &Contacts::is_sorted)
      .def_readwrite("instance_name", &Contacts::instance_name)

      .def("set_luni",             &Contacts::set_luni, py::arg("luni"))
      .def("sort_interactions",    &Contacts::sort_interactions)
      .def("sort_if_not_sorted",   &Contacts::sort_if_not_sorted)
      .def_static("prepare_input", &Contacts::prepare_input, py::arg("lhs"), py::arg("rhs"))

      .def("add", py::overload_cast<const Contacts &>(&Contacts::add), py::arg("contacts"))
      .def("add", py::overload_cast<const Contact &> (&Contacts::add), py::arg("interaction"))
      .def("add", py::overload_cast<EntityID, EntityID, float, InteractionType>(&Contacts::add))
      .def("add", py::overload_cast<const std::vector<Contact> &>(&Contacts::add), py::arg("interactions"))

      .def("add_many", py::overload_cast<const NSResults&, const std::vector<EntityID>&,                               InteractionType>(&Contacts::add_many))
      .def("add_many", py::overload_cast<const NSResults&, const std::vector<EntityID>&, const std::vector<EntityID>&, InteractionType>(&Contacts::add_many))

      .def("size",         &Contacts::size)
      .def("make_generic", &Contacts::make_generic)
      .def("get_luni",     &Contacts::get_luni, py::return_value_policy::reference)
      .def("print",        &Contacts::print_interactions)

      .def("__eq__",   &Contacts::operator==, py::is_operator())
      .def("__ne__",   &Contacts::operator!=, py::is_operator())
      .def("__and__",  [](Contacts &self, Contacts &other) {return self &  other;}, py::is_operator())
      .def("__or__",   [](Contacts &self, Contacts &other) {return self |  other;}, py::is_operator())
      .def("__sub__",  [](Contacts &self, Contacts &other) {return self -  other;}, py::is_operator())
      .def("__xor__",  [](Contacts &self, Contacts &other) {return self ^  other;}, py::is_operator())
      .def("__iand__", [](Contacts &self, Contacts &other) {return self &= other;}, py::is_operator())
      .def("__ior__",  [](Contacts &self, Contacts &other) {return self |= other;}, py::is_operator())
      .def("__isub__", [](Contacts &self, Contacts &other) {return self -= other;}, py::is_operator())
      .def("__ixor__", [](Contacts &self, Contacts &other) {return self ^= other;}, py::is_operator())

       // FIX: add support for slicing
      .def("__getitem__", [](const Contacts &self, int index) -> Contact {
        if (index < 0 || index >= self.size()) { throw py::index_error("Index out of range"); }
            return self.interactions[index];
       });
}
