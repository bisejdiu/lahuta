#ifndef LAHUTA_BINDINGS_INTERACTIONS_HPP
#define LAHUTA_BINDINGS_INTERACTIONS_HPP

#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>

#include "entities/interaction_types.hpp"

// clang-format off
namespace py = pybind11;
namespace lahuta::bindings {

inline void add_python_interactions(InteractionTypeSet& out, py::handle obj) {
  if (!obj || obj.is_none()) return;

  if (py::isinstance<InteractionTypeSet>(obj)) {
    out |= obj.cast<InteractionTypeSet>();
    return;
  }

  if (py::isinstance<InteractionType>(obj)) {
    out |= obj.cast<InteractionType>();
    return;
  }

  if (PySequence_Check(obj.ptr())) {
    py::sequence seq = py::reinterpret_borrow<py::sequence>(obj);
    for (auto item : seq) {
      add_python_interactions(out, item);
    }
    return;
  }

  throw py::type_error("Expected InteractionType, InteractionTypeSet, or iterable of interaction types");
}

} // namespace lahuta::bindings

#endif // LAHUTA_BINDINGS_INTERACTIONS_HPP
