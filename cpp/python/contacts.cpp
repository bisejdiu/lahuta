#include "contacts.hpp"
#include "lahuta.hpp"
#include "nn.hpp"

using namespace lahuta;

void bind_contacts(py::module &m) {

  py::enum_<InteractionType>(m, "InteractionType")
      .value("None", InteractionType::None)
      .value("Any", InteractionType::Any)
      .value("Hydrophobic", InteractionType::Hydrophobic)
      .value("HydrogenBond", InteractionType::HydrogenBond)
      .value("WeakHydrogenBond", InteractionType::WeakHydrogenBond)
      .value("Ionic", InteractionType::Ionic)
      .export_values();

  py::enum_<lahuta::EntityType>(m, "EntityType")
      .value("Atom", lahuta::EntityType::Atom)
      .value("Ring", lahuta::EntityType::Ring)
      .value("Group", lahuta::EntityType::Group)
      .export_values();

  py::class_<EntityID>(m, "EntityID") //
      .def(py::init<uint64_t>())
      .def("__int__", [](const EntityID &self) { return static_cast<uint64_t>(self); });

  py::class_<Contact>(m, "Contact")
      .def(py::init<EntityID, EntityID, float>(), py::arg("index1"), py::arg("index2"), py::arg("distance"))
      .def(
          py::init<EntityID, EntityID, float, InteractionType>(),
          py::arg("index1"),
          py::arg("index2"),
          py::arg("distance"),
          py::arg("type"))
      .def_readwrite("index1", &Contact::entity1)
      .def_readwrite("index2", &Contact::entity2)
      .def_readwrite("distance", &Contact::distance)
      .def_readwrite("type", &Contact::type)

      // Overloaded operators
      .def("__eq__", &Contact::operator==, py::is_operator())
      .def("__ne__", &Contact::operator!=, py::is_operator())
      .def("__lt__", &Contact::operator<, py::is_operator())
      .def("__gt__", &Contact::operator>, py::is_operator())
      .def("__le__", &Contact::operator<=, py::is_operator())
      .def("__ge__", &Contact::operator>=, py::is_operator())

      .def("__repr__", [](const Contact &self) {
        return "<Contact index1=" + std::to_string(self.entity1) //
               + ", index2=" + std::to_string(self.entity2)      //
               + ", distance=" + std::to_string(self.distance)
               + ", type=" + std::to_string(static_cast<int>(self.type)) + ">";
      });

  py::class_<Contacts>(m, "Contacts")
      // Constructors
      .def(py::init<>())
      .def(py::init<const class Luni *>(), py::arg("luni"))

      // Member variables
      .def_readwrite("interactions", &Contacts::interactions)
      .def_readwrite("is_sorted", &Contacts::is_sorted)
      .def_readwrite("instance_name", &Contacts::instance_name)

      .def("set_luni", &Contacts::set_luni, py::arg("luni"))
      .def("sort_interactions", &Contacts::sort_interactions)
      .def("sort_if_not_sorted", &Contacts::sort_if_not_sorted)
      .def_static("prepare_input", &Contacts::prepare_input, py::arg("lhs"), py::arg("rhs"))
      .def("add", py::overload_cast<const Contacts &>(&Contacts::add), py::arg("contacts"))
      .def("add", py::overload_cast<const Contact &>(&Contacts::add), py::arg("interaction"))
      .def("add", py::overload_cast<EntityID, EntityID, float, InteractionType>(&Contacts::add),
          py::arg("e1"),
          py::arg("e2"),
          py::arg("d"),
          py::arg("t"))
      .def("add", py::overload_cast<const std::vector<Contact> &>(&Contacts::add), py::arg("interactions"))
      .def("add_many",
          (void(Contacts::*)(
              const NSResults &,
              const std::vector<EntityID> &,
              const std::vector<EntityID> &,
              InteractionType))
              & Contacts::add_many,
          py::arg("neighbors"),
          py::arg("e1"),
          py::arg("e2"),
          py::arg("type") = InteractionType::Any)
      .def("add_many",
          (void(Contacts::*)(const NSResults &, const std::vector<EntityID> &, InteractionType))
              & Contacts::add_many,
          py::arg("neighbors"),
          py::arg("entities"),
          py::arg("type") = InteractionType::Any)
      .def("make_generic", &Contacts::make_generic)
      .def("size", &Contacts::size)
      .def("get_luni", &Contacts::get_luni, py::return_value_policy::reference)
      .def("print_interactions", &Contacts::print_interactions)

      // Overloaded operators
      .def("__eq__", &Contacts::operator==, py::is_operator())
      .def("__ne__", &Contacts::operator!=, py::is_operator())
      .def("__and__",  [](Contacts &self, Contacts &other) {return self &  other;}, py::is_operator())
      .def("__or__",   [](Contacts &self, Contacts &other) {return self |  other;}, py::is_operator())
      .def("__sub__",  [](Contacts &self, Contacts &other) {return self -  other;}, py::is_operator())
      .def("__xor__",  [](Contacts &self, Contacts &other) {return self ^  other;}, py::is_operator())
      .def("__iand__", [](Contacts &self, Contacts &other) {return self &= other;}, py::is_operator())
      .def("__ior__",  [](Contacts &self, Contacts &other) {return self |= other;}, py::is_operator())
      .def("__isub__", [](Contacts &self, Contacts &other) {return self -= other;}, py::is_operator())
      .def("__ixor__", [](Contacts &self, Contacts &other) {return self ^= other;}, py::is_operator())

      // Indexing support
      .def("__getitem__", [](const Contacts &self, int index) -> Contact {
            if (index < 0 || index >= self.size()) {
              throw py::index_error("Index out of range");
            }
            return self.interactions[index];
          },
          py::arg("index")
        );

  py::class_<InteractionOptions>(m, "InteractionOptions") //
      .def(py::init<float>(), py::arg("distance_cutoff"));

  py::class_<Interactions>(m, "Interactions")

      .def(py::init<Luni &, InteractionOptions>(), py::arg("luni"), py::arg("opts"))
      .def("hbond", &Interactions::hbond)
      .def("weak_hbond", &Interactions::weak_hbond)
      .def("hydrophobic", &Interactions::hydrophobic)
      .def("halogen", &Interactions::halogen)
      .def("ionic", &Interactions::ionic)
      .def("metalic", &Interactions::metalic)
      .def("cationpi", &Interactions::cationpi)
      .def("pistacking", &Interactions::pistacking);
}
