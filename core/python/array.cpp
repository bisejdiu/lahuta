#include "array.hpp"
#include "struct_unit.hpp"

using namespace lahuta;

py::array_t<float> coordinates(const RDGeom::POINT3D_VECT &coords) {
  if (coords.empty() || coords[0].dimension() != 3) {
    throw std::runtime_error("Invalid input: expected non-empty vector of 3D coordinates");
  }

  size_t n_atoms = coords.size();
  size_t n_dims = coords[0].dimension();

  auto result = py::array_t<float>({n_atoms, n_dims});
  auto buf = result.request();
  float *ptr = static_cast<float *>(buf.ptr);

  for (size_t i = 0; i < n_atoms; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      ptr[i * 3 + j] = coords[i][j];
    }
  }

  return result;
}

py::array point3d_to_pyarray(const RDGeom::Point3D &vec) {
  auto result = py::array_t<float>(3);
  auto buf = result.request();
  float *ptr = static_cast<float *>(buf.ptr);
  ptr[0] = vec.x;
  ptr[1] = vec.y;
  ptr[2] = vec.z;
  return result;
}

py::array string_array(const std::vector<std::string> &data) {
  size_t n = data.size();
  // creates an array of object dtype
  auto result = py::array_t<PyObject *>(n);

  auto buf = result.request();
  auto **ptr = static_cast<PyObject **>(buf.ptr);

  for (size_t i = 0; i < n; ++i) {
    ptr[i] = PyUnicode_FromString(data[i].c_str());
  }

  return result;
}

py::array int_array(const std::vector<int> &data) {
  size_t n = data.size();
  auto result = py::array_t<int>(n);
  auto buf = result.request();
  int *ptr = static_cast<int *>(buf.ptr);

  for (size_t i = 0; i < n; ++i) {
    ptr[i] = data[i];
  }

  return result;
}

py::array float_array(const std::vector<float> &data) {
  size_t n = data.size();
  auto result = py::array_t<float>(n);
  auto buf = result.request();
  float *ptr = static_cast<float *>(buf.ptr);

  for (size_t i = 0; i < n; ++i) {
    ptr[i] = data[i];
  }

  return result;
}

py::tuple factorize_residues(
    const std::vector<std::string> &resnames, const std::vector<int> &resids,
    const std::vector<std::string> &chains) {

  auto result = Factorizer::factorize({resnames, resids, chains});

  return py::make_tuple(
      py::array_t<int>(result.indices.size(), result.indices.data()),
      py::cast(result.resnames),
      py::cast(result.resids),
      py::cast(result.chainlabels));
}


void bind_common(py::module &m) {
  m.def("factorize_residues", &factorize_residues,
        "Factorize residue information into unique combinations");
}
