#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "lahuta.hpp"
#include "properties/analyzer.hpp"
#include "properties/query.hpp"
#include "properties/types.hpp"

// clang-format off
// using namespace lahuta;
namespace py = pybind11;

namespace {
using lahuta::PropertyKey;
}

namespace lahuta::bindings {

void bind_properties(py::module_ &m) {

  py::enum_<PropertyKey> PropertyKey_(m, "PropertyKey");

  py::class_<PropertyQuery<Luni>>     PropertyQuery_    (m, "PropertyQueryLuni");
  py::class_<PropertyAnalyzer<Luni>>  PropertyAnalyzer_ (m, "PropertyAnalyzerLuni");
  py::class_<IdentityAnalyzer<Luni>>  IdentityAnalyzer_ (m, "IdentityAnalyzerLuni");


  PropertyKey_
      .value("Names",     PropertyKey::Names)
      .value("Indices",   PropertyKey::Indices)
      .value("Elements",  PropertyKey::Elements)
      .value("Positions", PropertyKey::Positions);

  PropertyQuery_
      .def(py::init<>())
      .def("select", [](PropertyQuery<Luni> &self, PropertyKey key) -> PropertyQuery<Luni> & { return self.select(key); })
      .def("select", [](PropertyQuery<Luni> &self, const std::vector<PropertyKey> &keys) -> PropertyQuery<Luni> & {
            for (auto &key : keys) {
              self.select(key);
            }
            return self;
          })
      .def("properties", &PropertyQuery<Luni>::properties, py::return_value_policy::reference_internal);

  PropertyAnalyzer_
      .def(py::init<PropertyQuery<Luni>>());

  IdentityAnalyzer_
      .def(py::init<>());
}

} // namespace lahuta::bindings
