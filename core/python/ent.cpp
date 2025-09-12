#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "entities/entity_id.hpp"
#include "entities/find_contacts.hpp"
#include "entities/interaction_types.hpp"
#include "entities/records.hpp"
#include "topology.hpp"

namespace py = pybind11;
using namespace lahuta;

struct _PyEmptyParams final {}; // so we can create a ContactContext without requiring specific params

// clang-format off
namespace {

template <typename Rec> static inline auto make_predicate(py::function f) {
  return [f](const Rec &rec) { return f(rec).template cast<bool>(); };
}

static inline auto make_tester(py::function f) {
  return [f](uint32_t i, uint32_t j, float d2, const ContactContext &) -> InteractionType {
    if (!f || f.is_none()) return InteractionType::Generic;

    py::object obj = f(i, j, d2);
    if (py::isinstance<py::bool_>(obj)) {
      return obj.cast<bool>() ? InteractionType::Generic : InteractionType::None;
    }

    return obj.cast<InteractionType>();
  };
}

template <typename RecA, typename RecB>
static inline ContactSet dispatch_dual(
  const Topology &top,
  py::function pred_a, py::function pred_b,
  const search::SearchOptions &opts, py::function tester_py) {
  auto pa     = make_predicate<RecA>(pred_a);
  auto pb     = make_predicate<RecB>(pred_b);
  auto tester = make_tester(tester_py);
  _PyEmptyParams params;
  ContactContext ctx(top, params);
  return find_contacts(ctx, pa, pb, opts, tester);
}

template <typename Rec>
static inline ContactSet dispatch_self(
  const Topology &top,
  py::function pred,
  const search::SearchOptions &opts, py::function tester_py) {
  auto p      = make_predicate<Rec>(pred);
  auto tester = make_tester(tester_py);
  _PyEmptyParams params;
  ContactContext ctx(top, params);
  return find_contacts(ctx, p, opts, tester);
}

static inline ContactSet route_dual(
  const Topology &top,
  Kind kind_a, py::function pred_a,
  Kind kind_b, py::function pred_b,
  const search::SearchOptions &opts, py::function tester_py) {
  switch (kind_a) {
    case Kind::Atom:
      switch (kind_b) {
        case Kind::Atom:  return dispatch_dual<AtomRec, AtomRec> (top, pred_a, pred_b, opts, tester_py);
        case Kind::Ring:  return dispatch_dual<AtomRec, RingRec> (top, pred_a, pred_b, opts, tester_py);
        case Kind::Group: return dispatch_dual<AtomRec, GroupRec>(top, pred_a, pred_b, opts, tester_py);
      }
      break;
    case Kind::Ring:
      switch (kind_b) {
        case Kind::Atom:  return dispatch_dual<RingRec, AtomRec> (top, pred_a, pred_b, opts, tester_py);
        case Kind::Ring:  return dispatch_dual<RingRec, RingRec> (top, pred_a, pred_b, opts, tester_py);
        case Kind::Group: return dispatch_dual<RingRec, GroupRec>(top, pred_a, pred_b, opts, tester_py);
      }
      break;
    case Kind::Group:
      switch (kind_b) {
        case Kind::Atom:  return dispatch_dual<GroupRec, AtomRec> (top, pred_a, pred_b, opts, tester_py);
        case Kind::Ring:  return dispatch_dual<GroupRec, RingRec> (top, pred_a, pred_b, opts, tester_py);
        case Kind::Group: return dispatch_dual<GroupRec, GroupRec>(top, pred_a, pred_b, opts, tester_py);
      }
      break;
  }
  throw std::invalid_argument("Invalid Kind combination for find_contacts");
}

static inline ContactSet route_self(const Topology &top, Kind kind, py::function pred, const search::SearchOptions &opts, py::function tester_py) {
  switch (kind) {
    case Kind::Atom:  return dispatch_self<AtomRec> (top, pred, opts, tester_py);
    case Kind::Ring:  return dispatch_self<RingRec> (top, pred, opts, tester_py);
    case Kind::Group: return dispatch_self<GroupRec>(top, pred, opts, tester_py);
  }
  throw std::invalid_argument("Invalid Kind for self find_contacts");
}

} // namespace

void bind_entities(py::module &m) {
  py::enum_<FeatureGroup>(m, "FeatureGroup")
    .value("None",            FeatureGroup::None)
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

  py::class_<search::SearchOptions>(m, "SearchOptions")
    .def(py::init<>())
    .def_readwrite("distance_max",        &search::SearchOptions::distance_max)
    .def_readwrite("hit_reserve_factor",  &search::SearchOptions::hit_reserve_factor)
    .def_readwrite("sel_reserve_factor_a",&search::SearchOptions::sel_reserve_factor_a)
    .def_readwrite("sel_reserve_factor_b",&search::SearchOptions::sel_reserve_factor_b)
    .def_static("for_category", [](Category c) { return search::make_search_opts(c); }, py::arg("category"));

  // High-level contact discovery
  m.def("find_contacts",
        [](const Topology &top,
           Kind kind_a, py::function pred_a,
           Kind kind_b, py::function pred_b,
           const search::SearchOptions &opts,
           py::function tester) {
            return route_dual(top, kind_a, pred_a, kind_b, pred_b, opts, tester);
        },
        py::arg("topology"),
        py::arg("kind_a"), py::arg("pred_a"),
        py::arg("kind_b"), py::arg("pred_b"),
        py::arg("opts") = search::SearchOptions{},
        py::arg("tester") = py::none());

  m.def("find_contacts",
        [](const Topology &top,
           Kind kind, py::function pred,
           const search::SearchOptions &opts,
           py::function tester) {
            return route_self(top, kind, pred, opts, tester);
        },
        py::arg("topology"),
        py::arg("kind"), py::arg("pred"),
        py::arg("opts") = search::SearchOptions{},
        py::arg("tester") = py::none());
}
