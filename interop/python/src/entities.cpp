/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::array<const char*, 3> parts{"besian", "sejdiu", "@gmail.com"}; std::string s;
 *   for (auto p : parts) s.append(p, std::strlen(p));
 *   return s;
 * }();
 *
 */

#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/typing.h>

#include "compute/topology_snapshot.hpp"
#include "entities/contact.hpp"
#include "entities/context.hpp"
#include "entities/entity_id.hpp"
#include "entities/find_contacts.hpp"
#include "entities/interaction_types.hpp"
#include "entities/records.hpp"
#include "entities/resolver.hpp"
#include "entities/search/config.hpp"
#include "topology.hpp"

namespace py = pybind11;

struct _PyEmptyParams final {}; // To create a ContactContext without requiring specific params

namespace {
namespace L = lahuta;
namespace C = lahuta::compute;
using lahuta::search::SearchOptions;

template <typename Rec>
static inline std::function<bool(const Rec &)> make_predicate(const py::function &f) {
  // If no predicate is provided, select-all
  if (!f || f.is_none()) {
    return std::function<bool(const Rec &)>{[](const Rec &) { return true; }};
  }
  // Otherwise, call back into Python, acquire GIL
  return std::function<bool(const Rec &)>{[f](const Rec &rec) {
    py::gil_scoped_acquire ag;
    return f(rec).template cast<bool>();
  }};
}

static inline auto make_tester(const py::function &f) {
  return [f](uint32_t i, uint32_t j, float d2, const L::ContactContext &) -> L::InteractionType {
    if (!f || f.is_none()) return L::InteractionType::Generic;

    py::gil_scoped_acquire ag;
    py::object obj = f(i, j, d2);
    if (py::isinstance<py::bool_>(obj)) {
      return obj.cast<bool>() ? L::InteractionType::Generic : L::InteractionType::None;
    }
    return obj.cast<L::InteractionType>();
  };
}

template <typename RecA, typename RecB>
static inline L::ContactSet dispatch_dual(const lahuta::Topology &top, py::function pred_a,
                                          py::function pred_b, const lahuta::search::SearchOptions &opts,
                                          py::function tester_py) {
  auto pa     = make_predicate<RecA>(pred_a);
  auto pb     = make_predicate<RecB>(pred_b);
  auto tester = make_tester(tester_py);
  _PyEmptyParams params;
  auto tf = C::snapshot_of(top, top.conformer());
  L::ContactContext ctx(tf, params);
  return L::find_contacts(ctx, pa, pb, opts, tester);
}

template <typename Rec>
static inline L::ContactSet dispatch_self(const lahuta::Topology &top, py::function pred,
                                          const lahuta::search::SearchOptions &opts, py::function tester_py) {
  auto p      = make_predicate<Rec>(pred);
  auto tester = make_tester(tester_py);
  _PyEmptyParams params;
  auto tf = C::snapshot_of(top, top.conformer());
  L::ContactContext ctx(tf, params);
  return L::find_contacts(ctx, p, opts, tester);
}

static inline L::ContactSet route_dual(const lahuta::Topology &top, L::Kind kind_a, py::function pred_a,
                                       L::Kind kind_b, py::function pred_b,
                                       const lahuta::search::SearchOptions &opts, py::function tester_py) {
  // If there are no Python callbacks, release the GIL during the heavy C++ work
  const bool no_py_calls = (!pred_a || pred_a.is_none()) && (!pred_b || pred_b.is_none()) &&
                           (!tester_py || tester_py.is_none());
  std::unique_ptr<py::gil_scoped_release> gil_rel;
  if (no_py_calls) gil_rel = std::make_unique<py::gil_scoped_release>();
  switch (kind_a) {
    case L::Kind::Atom:
      switch (kind_b) {
        case L::Kind::Atom:
          return dispatch_dual<L::AtomRec, L::AtomRec>(top, pred_a, pred_b, opts, tester_py);
        case L::Kind::Ring:
          return dispatch_dual<L::AtomRec, L::RingRec>(top, pred_a, pred_b, opts, tester_py);
        case L::Kind::Group:
          return dispatch_dual<L::AtomRec, L::GroupRec>(top, pred_a, pred_b, opts, tester_py);
      }
      break;
    case L::Kind::Ring:
      switch (kind_b) {
        case L::Kind::Atom:
          return dispatch_dual<L::RingRec, L::AtomRec>(top, pred_a, pred_b, opts, tester_py);
        case L::Kind::Ring:
          return dispatch_dual<L::RingRec, L::RingRec>(top, pred_a, pred_b, opts, tester_py);
        case L::Kind::Group:
          return dispatch_dual<L::RingRec, L::GroupRec>(top, pred_a, pred_b, opts, tester_py);
      }
      break;
    case L::Kind::Group:
      switch (kind_b) {
        case L::Kind::Atom:
          return dispatch_dual<L::GroupRec, L::AtomRec>(top, pred_a, pred_b, opts, tester_py);
        case L::Kind::Ring:
          return dispatch_dual<L::GroupRec, L::RingRec>(top, pred_a, pred_b, opts, tester_py);
        case L::Kind::Group:
          return dispatch_dual<L::GroupRec, L::GroupRec>(top, pred_a, pred_b, opts, tester_py);
      }
      break;
  }
  throw std::invalid_argument("Invalid Kind combination for find_contacts");
}

static inline L::ContactSet route_self(const lahuta::Topology &top, L::Kind kind, py::function pred,
                                       const lahuta::search::SearchOptions &opts, py::function tester_py) {
  const bool no_py_calls = (!pred || pred.is_none()) && (!tester_py || tester_py.is_none());
  std::unique_ptr<py::gil_scoped_release> gil_rel;
  if (no_py_calls) gil_rel = std::make_unique<py::gil_scoped_release>();
  switch (kind) {
    case L::Kind::Atom:
      return dispatch_self<L::AtomRec>(top, pred, opts, tester_py);
    case L::Kind::Ring:
      return dispatch_self<L::RingRec>(top, pred, opts, tester_py);
    case L::Kind::Group:
      return dispatch_self<L::GroupRec>(top, pred, opts, tester_py);
  }
  throw std::invalid_argument("Invalid Kind for self find_contacts");
}

} // namespace

namespace lahuta::bindings {

void bind_entities(py::module &m) {

  py::class_<SearchOptions>(m, "SearchOptions")
      .def(py::init<>())
      .def_readwrite("distance_max", &SearchOptions::distance_max);

  // Dual-entity overload without tester
  m.def(
      "find_contacts",
      [](const lahuta::Topology &top,
         L::Kind kind_a,
         py::function pred_a,
         L::Kind kind_b,
         py::function pred_b,
         const SearchOptions &opts) {
        return route_dual(top, kind_a, pred_a, kind_b, pred_b, opts, py::function{});
      },
      py::arg("topology"),
      py::arg("kind_a"),
      py::arg("pred_a"),
      py::arg("kind_b"),
      py::arg("pred_b"),
      py::arg_v("opts", SearchOptions{}, "SearchOptions()"));

  // Dual-entity overload with tester (tester before opts to avoid non-default after default in stubs)
  m.def(
      "find_contacts",
      [](const lahuta::Topology &top,
         L::Kind kind_a,
         py::function pred_a,
         L::Kind kind_b,
         py::function pred_b,
         py::function tester,
         const SearchOptions &opts) { return route_dual(top, kind_a, pred_a, kind_b, pred_b, opts, tester); },
      py::arg("topology"),
      py::arg("kind_a"),
      py::arg("pred_a"),
      py::arg("kind_b"),
      py::arg("pred_b"),
      py::arg("tester"),
      py::arg_v("opts", SearchOptions{}, "SearchOptions()"));

  // Self-entity overload without tester
  m.def(
      "find_contacts",
      [](const lahuta::Topology &top, L::Kind kind, py::function pred, const SearchOptions &opts) {
        return route_self(top, kind, pred, opts, py::function{});
      },
      py::arg("topology"),
      py::arg("kind"),
      py::arg("pred"),
      py::arg_v("opts", SearchOptions{}, "SearchOptions()"));

  // Self-entity overload with tester (tester before opts to avoid non-default after default in stubs)
  m.def(
      "find_contacts",
      [](const lahuta::Topology &top,
         L::Kind kind,
         py::function pred,
         py::function tester,
         const SearchOptions &opts) { return route_self(top, kind, pred, opts, tester); },
      py::arg("topology"),
      py::arg("kind"),
      py::arg("pred"),
      py::arg("tester"),
      py::arg_v("opts", SearchOptions{}, "SearchOptions()"));

  // EntityResolver: thin facade over Topology::resolve<K>, providing ergonomic APIs to Python.
  py::class_<L::EntityResolver>(m, "EntityResolver", R"doc(
Resolve EntityID to concrete records (AtomRec, RingRec, GroupRec) without copying.

Notes
- The resolver holds a reference to the provided Topology; the topology is kept
  alive for as long as the resolver exists (via keep_alive).
- Returned records are references to internal storage; do not store them beyond
  the lifetime of the topology.
)doc")
      .def(py::init<const lahuta::Topology &>(),
           py::arg("topology"),
           py::keep_alive<1, 2>(),
           R"doc(Create a resolver bound to a specific topology.)doc")
      .def_property_readonly("topology",
                             &L::EntityResolver::topology,
                             py::return_value_policy::reference_internal,
                             R"doc(Underlying topology (borrowed reference).)doc")

      .def(
          "resolve",
          [](L::EntityResolver &self, const L::EntityID &id) -> py::object {
            const auto &top = self.topology();
            switch (id.kind()) {
              case L::Kind::Atom:
                return py::cast(&top.resolve<L::Kind::Atom>(id),
                                py::return_value_policy::reference,
                                py::cast(&top, py::return_value_policy::reference));
              case L::Kind::Ring:
                return py::cast(&top.resolve<L::Kind::Ring>(id),
                                py::return_value_policy::reference,
                                py::cast(&top, py::return_value_policy::reference));
              case L::Kind::Group:
                return py::cast(&top.resolve<L::Kind::Group>(id),
                                py::return_value_policy::reference,
                                py::cast(&top, py::return_value_policy::reference));
              default:
                throw std::invalid_argument("Invalid EntityID kind");
            }
          },
          py::arg("id"),
          R"doc(Resolve an EntityID to a concrete record. Returns AtomRec, RingRec, or GroupRec.)doc")

      //
      // Return type is a union of all 9 precise tuple combinations
      // Obviously, I'm not happy with this, but it's the best *I* know how to do. - Besian, August 2025
      //
      .def(
          "resolve_contact",
          [](L::EntityResolver &self,
             const L::Contact &c) -> py::typing::Union<py::typing::Tuple<L::AtomRec, L::AtomRec>,
                                                       py::typing::Tuple<L::AtomRec, L::RingRec>,
                                                       py::typing::Tuple<L::AtomRec, L::GroupRec>,
                                                       py::typing::Tuple<L::RingRec, L::AtomRec>,
                                                       py::typing::Tuple<L::RingRec, L::RingRec>,
                                                       py::typing::Tuple<L::RingRec, L::GroupRec>,
                                                       py::typing::Tuple<L::GroupRec, L::AtomRec>,
                                                       py::typing::Tuple<L::GroupRec, L::RingRec>,
                                                       py::typing::Tuple<L::GroupRec, L::GroupRec>> {
            const auto &top = self.topology();
            py::object a, b;
            switch (c.lhs.kind()) {
              case L::Kind::Atom:
                a = py::cast(&top.resolve<L::Kind::Atom>(c.lhs),
                             py::return_value_policy::reference,
                             py::cast(&top, py::return_value_policy::reference));
                break;
              case L::Kind::Ring:
                a = py::cast(&top.resolve<L::Kind::Ring>(c.lhs),
                             py::return_value_policy::reference,
                             py::cast(&top, py::return_value_policy::reference));
                break;
              case L::Kind::Group:
                a = py::cast(&top.resolve<L::Kind::Group>(c.lhs),
                             py::return_value_policy::reference,
                             py::cast(&top, py::return_value_policy::reference));
                break;
              default:
                throw std::invalid_argument("Invalid lhs kind");
            }
            switch (c.rhs.kind()) {
              case L::Kind::Atom:
                b = py::cast(&top.resolve<L::Kind::Atom>(c.rhs),
                             py::return_value_policy::reference,
                             py::cast(&top, py::return_value_policy::reference));
                break;
              case L::Kind::Ring:
                b = py::cast(&top.resolve<L::Kind::Ring>(c.rhs),
                             py::return_value_policy::reference,
                             py::cast(&top, py::return_value_policy::reference));
                break;
              case L::Kind::Group:
                b = py::cast(&top.resolve<L::Kind::Group>(c.rhs),
                             py::return_value_policy::reference,
                             py::cast(&top, py::return_value_policy::reference));
                break;
              default:
                throw std::invalid_argument("Invalid rhs kind");
            }
            return py::make_tuple(std::move(a), std::move(b));
          },
          py::arg("contact"),
          R"doc(Resolve both sides of a contact. Returns a tuple (left_record, right_record).)doc")

      // Return type is list[utple[AtomRec|RingRec|GroupRec, AtomRec|RingRec|GroupRec]]
      .def(
          "resolve_all",
          [](L::EntityResolver &self, const L::ContactSet &set)
              -> py::typing::List<py::typing::Tuple<py::typing::Union<L::AtomRec, L::RingRec, L::GroupRec>,
                                                    py::typing::Union<L::AtomRec, L::RingRec, L::GroupRec>>> {
            const auto &top = self.topology();
            py::list out;
            for (const auto &c : set.data()) {
              py::object a, b;
              switch (c.lhs.kind()) {
                case L::Kind::Atom:
                  a = py::cast(&top.resolve<L::Kind::Atom>(c.lhs),
                               py::return_value_policy::reference,
                               py::cast(&top, py::return_value_policy::reference));
                  break;
                case L::Kind::Ring:
                  a = py::cast(&top.resolve<L::Kind::Ring>(c.lhs),
                               py::return_value_policy::reference,
                               py::cast(&top, py::return_value_policy::reference));
                  break;
                case L::Kind::Group:
                  a = py::cast(&top.resolve<L::Kind::Group>(c.lhs),
                               py::return_value_policy::reference,
                               py::cast(&top, py::return_value_policy::reference));
                  break;
                default:
                  throw std::invalid_argument("Invalid lhs kind");
              }
              switch (c.rhs.kind()) {
                case L::Kind::Atom:
                  b = py::cast(&top.resolve<L::Kind::Atom>(c.rhs),
                               py::return_value_policy::reference,
                               py::cast(&top, py::return_value_policy::reference));
                  break;
                case L::Kind::Ring:
                  b = py::cast(&top.resolve<L::Kind::Ring>(c.rhs),
                               py::return_value_policy::reference,
                               py::cast(&top, py::return_value_policy::reference));
                  break;
                case L::Kind::Group:
                  b = py::cast(&top.resolve<L::Kind::Group>(c.rhs),
                               py::return_value_policy::reference,
                               py::cast(&top, py::return_value_policy::reference));
                  break;
                default:
                  throw std::invalid_argument("Invalid rhs kind");
              }
              out.append(py::make_tuple(std::move(a), std::move(b)));
            }
            return out;
          },
          py::arg("contacts"),
          R"doc(Materialize a list of resolved record pairs for a ContactSet.)doc");
}
} // namespace lahuta::bindings
