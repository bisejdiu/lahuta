#ifndef LAHUTA_BINDINGS_HPP
#define LAHUTA_BINDINGS_HPP

#include <Python.h>
#include <pybind11/pybind11.h>
#if PY_VERSION_HEX < 0x030A0000
#error "Lahuta Python bindings require Python >= 3.10"
#endif

namespace py = pybind11;

namespace lahuta::bindings {

void bind_utilities(py::module_ &m);

void bind_luni(py::module_ &m);
void bind_topology(py::module_ &m);
void bind_records(py::module_ &m);
void bind_atom_types(py::module_ &m);
void bind_elements(py::module_ &m);
void bind_entity_id(py::module_ &m);
void bind_entities(py::module_ &m);

void bind_grid_search(py::module_ &m);
void bind_contacts(py::module_ &m);
void bind_distance(py::module_ &m);

void bind_foldseek(py::module_ &m);
void bind_align(py::module_ &m);
void bind_properties(py::module_ &m);
void bind_te(py::module_ &m);
void bind_rdkit(py::module_ &m);

void bind_logger(py::module_ &m);

void bind_pipeline_dynamic(py::module_ &m);
void bind_db(py::module_ &m);

} // namespace lahuta::bindings

#endif // LAHUTA_BINDINGS_HPP
