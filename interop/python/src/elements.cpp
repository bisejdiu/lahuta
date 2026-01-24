#include <string>

#include <pybind11/pybind11.h>

#include "bindings.hpp"
#include "chemistry/elements.hpp"

// clang-format off
namespace py = pybind11;
namespace {

using lahuta::Element;

inline Element element_from_atomic_number(unsigned int atomic_number) {
  constexpr unsigned int MinAtomic = 1;
  constexpr unsigned int MaxAtomic = static_cast<unsigned int>(Element::Og);

  if (atomic_number < MinAtomic) throw py::value_error("atomic_number must be >= 1");
  if (atomic_number > MaxAtomic) throw py::value_error("atomic_number exceeds supported range");

  return static_cast<Element>(atomic_number);
}

struct ElementWrapper {
  explicit ElementWrapper(Element element) : value(element) {}
  Element value;
};

ElementWrapper element_from_symbol(const std::string &symbol) {
  const auto element = lahuta::elements::find_element(symbol.c_str());
  const auto atomic_number = static_cast<unsigned int>(element);

  if (atomic_number == 0) throw py::value_error("Unknown element symbol: " + symbol);
  return ElementWrapper(element);
}

} // namespace

namespace lahuta::bindings {

void bind_elements(py::module_ &m) {
  auto element_cls = py::class_<ElementWrapper>(m, "Element", "Chemical element descriptor.");

  element_cls
      .def(py::init([](unsigned int atomic_number) {
             return ElementWrapper(element_from_atomic_number(atomic_number));
           }),
           py::arg("atomic_number"),
           "Create an element from its atomic number (Z).")
      .def_static(
          "from_symbol",
          [](const std::string &symbol) { return element_from_symbol(symbol); },
          py::arg("symbol"),
          "Create an element from its chemical symbol (case-sensitive).")
      .def_property_readonly(
          "atomic_number",
          [](const ElementWrapper &self) {
            return static_cast<unsigned int>(self.value);
          },
          "Atomic number (Z).")
      .def_property_readonly(
          "symbol",
          [](const ElementWrapper &self) {
            return std::string(::gemmi::Element(self.value).name());
          },
          "Chemical symbol for the element.")
      .def_property_readonly(
          "name",
          [](const ElementWrapper &self) {
            return std::string(lahuta::elements::element_name(self.value));
          },
          "IUPAC element name.")
      .def_property_readonly(
          "vdw_radius",
          [](const ElementWrapper &self) {
            return lahuta::elements::vdw_radius(self.value);
          },
          "van der Waals radius (A).")
      .def_property_readonly(
          "covalent_radius",
          [](const ElementWrapper &self) {
            return lahuta::elements::covalent_radius(self.value);
          },
          "Covalent radius (A).")
      .def("__repr__",
           [](const ElementWrapper &self) {
             const auto symbol = std::string(::gemmi::Element(self.value).name());
             return "Element(symbol='" + symbol + "', Z=" + std::to_string(static_cast<unsigned int>(self.value)) + ")";
           })
      .def("__str__", [](const ElementWrapper &self) {
        return std::string(lahuta::elements::element_name(self.value));
      });
  m.def(
      "vdw_radius",
      [](unsigned int atomic_number) {
        return lahuta::elements::vdw_radius(element_from_atomic_number(atomic_number));
      },
      py::arg("atomic_number"),
      "van der Waals radius (A) for the element with the given atomic number.");

  m.def(
      "covalent_radius",
      [](unsigned int atomic_number) {
        return lahuta::elements::covalent_radius(element_from_atomic_number(atomic_number));
      },
      py::arg("atomic_number"),
      "Covalent radius (A) for the element with the given atomic number.");
}

} // namespace lahuta::bindings
