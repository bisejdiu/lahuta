#include <pybind11/pybind11.h>

#include "bindings.hpp"

namespace py = pybind11;

PYBIND11_MODULE(mapping, m) {
    m.doc() = "Lahuta mapping submodule for sequence and structure alignment";

    py::options options;
    options.disable_enum_members_docstring();

    lahuta::bindings::bind_foldseek(m);
    lahuta::bindings::bind_align(m);
    lahuta::bindings::bind_te(m);
}
