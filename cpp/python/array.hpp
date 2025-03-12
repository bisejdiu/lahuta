#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <rdkit/Geometry/point.h>

namespace py = pybind11;

py::array_t<float> coordinates(const RDGeom::POINT3D_VECT &coords);
py::array point3d_to_pyarray(const RDGeom::Point3D &vec);
py::array string_array(const std::vector<std::string> &data);
py::array int_array(const std::vector<int> &data);
py::tuple factorize_residues(
    const std::vector<std::string> &resnames, const std::vector<int> &resids,
    const std::vector<std::string> &chains);

void bind_common(py::module &m);
