/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s;
 *   for (std::string_view part : {"besian", "sejdiu", "@gmail.com"})
 *     std::copy(part.begin(), part.end(), std::back_inserter(s));
 *   return s;
 * }();
 *
 */

#include <pybind11/pybind11.h>

#include "entities/entity_id.hpp"

namespace py = pybind11;

// clang-format off
namespace lahuta::bindings {

// Cache canonical Python Kind members to ensure identity (`is`) comparisons work.
static py::object PY_KIND_ATOM;
static py::object PY_KIND_RING;
static py::object PY_KIND_GROUP;

void bind_entity_id(py::module_ &m) {

  py::enum_<Kind>(m, "Kind", R"doc(
Entity category (Atom=0, Ring=1, Group=2).

Stored in the upper 8 bits of an EntityID and used for ordering.
)doc")
    .value("Atom",  Kind::Atom,  R"doc(Kind 0. Highest ordering priority among kinds.)doc")
    .value("Ring",  Kind::Ring,  R"doc(Kind 1. Intermediate ordering priority.)doc")
    .value("Group", Kind::Group, R"doc(Kind 2. Lowest ordering priority.)doc");

  // Capture canonical enum members so we can return them from C++
  try {
    py::object KindObj = m.attr("Kind");
    PY_KIND_ATOM  = KindObj.attr("Atom");
    PY_KIND_RING  = KindObj.attr("Ring");
    PY_KIND_GROUP = KindObj.attr("Group");
  } catch (...) {
    // If this fails, identity checks may not hold, but functionality remains.
  }

  py::class_<EntityID>(m, "EntityID", R"doc(
64-bit packed identifier: [ kind:8 | index:56 ].

- `kind` returns the canonical `Kind`
- `index` is a 32-bit view; `int(self)` yields the raw payload
- Ordering and hashing use the raw 64-bit value

Example
>>> from lahuta import Kind, EntityID
>>> eid = EntityID.make(Kind.Ring, 42)
>>> int(eid) == (1 << 56) | 42
True
>>> eid.kind is Kind.Ring
True
>>> eid.index
42
>>> str(eid)
'Ring#42'
)doc")
    .def(py::init<u64>(),
         py::arg("raw"),
         R"doc(
Construct from raw 64-bit payload `[ kind:8 | index:56 ]`.

No validation. `str(self)` raises `ValueError` if the kind is invalid.
)doc")

    .def_property_readonly(
        "kind",
        [](const EntityID &self) -> py::object {
          // Return canonical Python Kind members so `eid.kind is Kind.Ring` is True.
          switch (self.kind()) {
            case Kind::Atom:  return PY_KIND_ATOM;
            case Kind::Ring:  return PY_KIND_RING;
            case Kind::Group: return PY_KIND_GROUP;
            default:          return py::cast(self.kind());
          }
        },
        R"doc(
Kind encoded in the upper 8 bits; returns the canonical `Kind` (identity-safe).
)doc")

    .def_property_readonly(
        "index",
        &EntityID::index,
        R"doc(
Low 32-bit view of the 56-bit index. The full value is preserved in `int(self)`.
)doc")

    .def_static(
        "make",
        &EntityID::make,
        py::arg("kind"),
        py::arg("index"),
        R"doc(
Pack `(kind, index)` into an `EntityID`.

The raw 56-bit index is preserved; `.index` returns its low 32 bits.
)doc")

    .def("__int__",
         [](const EntityID &self) { return static_cast<std::uint64_t>(self.raw); },
         R"doc(
Exact 64-bit raw payload (`int`). Layout: `[ kind:8 | index:56 ]`.
)doc")

    .def("__str__",
         [](const EntityID &self) { return self.to_string(); },
         R"doc(
Human-readable form: `"Kind#index"`. Raises `ValueError` if the kind is invalid.
)doc")

    .def("__repr__",
         [](const EntityID &self) { return "EntityID(" + self.to_string() + ")"; },
         R"doc(
Debug representation: `"EntityID(Kind#index)"`.
)doc")

    .def(
        "__hash__",
        [](const EntityID &self) {
          return static_cast<std::size_t>(std::hash<u64>{}(self.raw));
        },
        R"doc(
Hash of the raw 64-bit payload (consistent with equality).
)doc")

    .def(
        "__eq__",
        [](const EntityID &self, py::object other) -> py::object {
          if (py::isinstance<EntityID>(other)) {
            const auto &o = other.cast<const EntityID &>();
            return py::bool_(self == o);
          }
          // Signal to Python that the comparison is not implemented for this type.
          return py::reinterpret_borrow<py::object>(Py_NotImplemented);
        },
        R"doc(
Equality by raw payload; returns `NotImplemented` for non-`EntityID` operands.
)doc")

    .def(
        "__ne__",
        [](const EntityID &self, py::object other) -> py::object {
          if (py::isinstance<EntityID>(other)) {
            const auto &o = other.cast<const EntityID &>();
            return py::bool_(self != o);
          }
          return py::reinterpret_borrow<py::object>(Py_NotImplemented);
        },
        R"doc(
Inequality by raw payload; returns `NotImplemented` for non-`EntityID` operands.
)doc")

    .def(
        "__lt__",
        [](const EntityID &self, py::object other) -> py::object {
          if (py::isinstance<EntityID>(other)) {
            const auto &o = other.cast<const EntityID &>();
            return py::bool_(self < o);
          }
          return py::reinterpret_borrow<py::object>(Py_NotImplemented);
        },
        R"doc(
Ordering by raw payload (equivalently by `(kind, index)`). Returns `NotImplemented` for other types.
)doc");
}
} // namespace lahuta::bindings
