#ifndef LAHUTA_PY_CONVERSIONS_HPP
#define LAHUTA_PY_CONVERSIONS_HPP

#include <cstring>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "distances.hpp"

// clang-format off
namespace py = pybind11;
namespace lahuta::numpy {

// Convert numpy (n,3) -> std::vector<std::vector<double>>
inline std::vector<std::vector<double>> array2vv_double(py::array_t<double, py::array::c_style | py::array::forcecast> arr) {
  if (arr.ndim()  != 2) throw std::invalid_argument("points must be a 2D array of shape (N, 3)");
  if (arr.shape(1) < 3) throw std::invalid_argument("points must have at least 3 columns (x,y,z)");

  const py::ssize_t N = arr.shape(0);
  const double *data = arr.data();

  std::vector<std::vector<double>> out;
  out.reserve(static_cast<size_t>(N));
  for (py::ssize_t i = 0; i < N; ++i) {
    const double *row = data + i * arr.shape(1);
    out.push_back({row[0], row[1], row[2]});
  }
  return out;
}

// Convert a Python sequence of sequences -> std::vector<std::vector<double>>
inline std::vector<std::vector<double>> seq2vv_double(const py::sequence &seq) {
  std::vector<std::vector<double>> out;
  out.reserve(seq.size());
  for (auto it : seq) {
    py::sequence row = py::cast<py::sequence>(it);
    if (row.size() < 3) throw std::invalid_argument("each point must have at least 3 coordinates (x,y,z)");
    out.push_back({py::cast<double>(row[0]), py::cast<double>(row[1]), py::cast<double>(row[2])});
  }
  return out;
}

// Copy ContiguousMatrix<T> into a new numpy array (rows, cols)
template <class T>
inline py::array_t<T> matrix_to_numpy(const lahuta::ContiguousMatrix<T> &M) {
  py::array_t<T> out({M.rows(), M.cols()});
  std::memcpy(out.mutable_data(), M.data(), sizeof(T) * static_cast<size_t>(M.rows()) * static_cast<size_t>(M.cols()));
  return out;
}

} // namespace lahuta::numpy

#endif // LAHUTA_PY_CONVERSIONS_HPP
