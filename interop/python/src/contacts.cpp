#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "analysis/contacts/provider.hpp"
#include "contacts/arpeggio/provider.hpp"
#include "contacts/engine.hpp"
#include "contacts/molstar/provider.hpp"

#include "entities/contact.hpp"
#include "entities/entity_id.hpp"
#include "entities/interaction_types.hpp"
#include "topology.hpp"

namespace py = pybind11;

// clang-format off
namespace lahuta::bindings {

void bind_contacts(py::module_ &m) {

  py::enum_<analysis::contacts::ContactProvider>(m, "ContactProvider")
    .value("MolStar",  analysis::contacts::ContactProvider::MolStar)
    .value("Arpeggio", analysis::contacts::ContactProvider::Arpeggio);

  py::enum_<Category>(m, "Category")
    .value("None_",                 Category::None)
    .value("Generic",               Category::Generic)
    .value("Hydrophobic",           Category::Hydrophobic)
    .value("Halogen",               Category::Halogen)
    .value("HydrogenBond",          Category::HydrogenBond)
    .value("WeakHydrogenBond",      Category::WeakHydrogenBond)
    .value("PolarHydrogenBond",     Category::PolarHydrogenBond)
    .value("WeakPolarHydrogenBond", Category::WeakPolarHydrogenBond)
    .value("Aromatic",              Category::Aromatic)
    .value("Ionic",                 Category::Ionic)
    .value("MetalCoordination",     Category::MetalCoordination)
    .value("CationPi",              Category::CationPi)
    .value("PiStacking",            Category::PiStacking)
    .value("Carbonyl",              Category::Carbonyl)
    .value("VanDerWaals",           Category::VanDerWaals)
    .value("DonorPi",               Category::DonorPi)
    .value("SulphurPi",             Category::SulphurPi)
    .value("CarbonPi",              Category::CarbonPi);

  py::enum_<Flavor>(m, "Flavor")
    .value("Default",  Flavor::Default)
    .value("Parallel", Flavor::Parallel)
    .value("TShape",   Flavor::TShape);

  // A concret type of molecular interaction with category and flavor
  py::class_<InteractionType>(m, "InteractionType")
    .def(py::init<Category, Flavor>(),
         py::arg_v("category", Category::None, "Category.None_"),
         py::arg_v("flavor",   Flavor::Default, "Flavor.Default"))
    .def_readonly_static("All",                   &InteractionType::All)
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

    .def("__str__",  [](const InteractionType& self) { return interaction_type_to_string(self); })
    .def("__repr__", [](const InteractionType& self) { return "InteractionType(" + interaction_type_to_string(self) + ")"; })
    .def(
      "__int__",
      [](const InteractionType& self) { return static_cast<unsigned int>(static_cast<std::uint32_t>(self)); },
      R"doc(
Return the 32-bit packed code for this InteractionType.

Layout: [ flavor:16 | category:16 ]. Both Category and Flavor are 16-bit enums.
)doc")

    .def("__hash__", [](const InteractionType& self) {
        return static_cast<std::size_t>(static_cast<std::uint32_t>(self));
      })
    .def("__eq__", [](const InteractionType& self, py::object other) -> py::object {
        if (py::isinstance<InteractionType>(other)) {
          auto o = other.cast<InteractionType>();
          return py::bool_(self == o);
        }
        return py::reinterpret_borrow<py::object>(Py_NotImplemented);
      })
    .def("__ne__", [](const InteractionType& self, py::object other) -> py::object {
        if (py::isinstance<InteractionType>(other)) {
          auto o = other.cast<InteractionType>();
          return py::bool_(self != o);
        }
        return py::reinterpret_borrow<py::object>(Py_NotImplemented);
      });

  // A contact between two entities with distance and interaction type
    py::class_<Contact>(m, "Contact")
      .def(py::init<EntityID, EntityID, float, InteractionType>(),
           py::arg("lhs"), py::arg("rhs"), py::arg("distance_sq"), py::arg("type"))
      .def_property(
        "lhs",
        [](const Contact &self) { return self.lhs; },
        [](Contact &self, const EntityID &value) { self.lhs = value; },
        "Left entity (EntityID)"
      )
      .def_property(
        "rhs",
        [](const Contact &self) { return self.rhs; },
        [](Contact &self, const EntityID &value) { self.rhs = value; },
        "Right entity (EntityID)"
      )
      .def_property("distance_sq",
        [](const Contact &self) { return self.distance; },
        [](Contact &self, float value) { self.distance = value; },
        "Distance squared between entities (Å^2)" )
      .def_property("type", 
        [](const Contact &self) { return self.type; },
        [](Contact &self, const InteractionType &value) { self.type = value; },
        "Interaction type (InteractionType)")

      .def("__hash__", [](const Contact& c) {
          size_t h1 = std::hash<uint64_t>{}(c.lhs.raw);
          size_t h2 = std::hash<uint64_t>{}(c.rhs.raw);
          size_t h3 = std::hash<uint32_t>{}(static_cast<uint32_t>(c.type));
          return (h1 ^ (h2 << 1)) ^ (h3 << 2);
      })
      .def("__eq__", [](const Contact &self, py::object other) -> bool {
        if (!py::isinstance<Contact>(other)) return false;
        return self == other.cast<Contact>();
      }, py::is_operator())
      .def("__ne__", [](const Contact &self, py::object other) -> bool {
        if (!py::isinstance<Contact>(other)) return true;
        return self != other.cast<Contact>();
      }, py::is_operator())
      .def("__lt__", &Contact::operator<,  py::is_operator())
      .def("__repr__", [](const Contact &self) {
        return py::str("Contact(lhs={}, rhs={}, type={}, distance_sq={})")
          .format(self.lhs.to_string(), self.rhs.to_string(), interaction_type_to_string(self.type), self.distance);
      })
      .def("__str__", [](const Contact &self) {
        return py::str("Contact(lhs={}, rhs={}, type={}, distance_sq={})")
          .format(self.lhs.to_string(), self.rhs.to_string(), interaction_type_to_string(self.type), self.distance);
      })
      .def("to_dict", [](const Contact &self) {
        py::dict d;
        d[py::str("lhs_kind")]    = self.lhs.kind();
        d[py::str("lhs_index")]   = self.lhs.index();
        d[py::str("rhs_kind")]    = self.rhs.kind();
        d[py::str("rhs_index")]   = self.rhs.index();
        d[py::str("distance_sq")] = self.distance;
        d[py::str("type")]        = self.type;
        d[py::str("category")]    = self.type.category;
        d[py::str("flavor")]      = self.type.flavor;
        return d;
      }, "Return a dict representation suitable for testing/serialization")
      ;

  // Ordered container for Contact obj with set operations and no duplicates
  py::class_<ContactSet>(m, "ContactSet")
    .def(py::init<>())
    .def("data",         [](ContactSet &self) { return self.data(); })
    .def("insert",       py::overload_cast<const ContactSet &>(&ContactSet::insert), py::arg("contacts"))
    .def("insert",       py::overload_cast<const Contact    &>(&ContactSet::insert), py::arg("contact" ))

    .def("set_union",         &ContactSet::set_union, py::arg("other"))
    .def("set_intersection",  &ContactSet::set_intersection, py::arg("other"))
    .def("set_difference",    &ContactSet::set_difference, py::arg("other"))
    .def("set_symmetric_difference", &ContactSet::set_symmetric_difference, py::arg("other"))

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
} // namespace lahuta::bindings
