#pragma once

#include <pybind11/pybind11.h>

namespace py = pybind11;

void xbind_atom_types(py::module &m);
