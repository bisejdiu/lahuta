#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "distances.hpp"
#include "distopia.h"
#include "py_conversions.hpp"

namespace py = pybind11;

// clang-format off
namespace {
template <typename T>
static py::array_t<T> distance(const std::vector<std::vector<T>> &points1, const std::vector<std::vector<T>> &points2) {
  const int na = static_cast<int>(points1.size());
  const int nb = static_cast<int>(points2.size());

  if (na == 0 || nb == 0) return py::array_t<T>();

  std::vector<T> a(3 * na);
  for (int i = 0; i < na; ++i) {
    a[3 * i + 0] = points1[i][0];
    a[3 * i + 1] = points1[i][1];
    a[3 * i + 2] = points1[i][2];
  }

  std::vector<T> b(3 * nb);
  for (int i = 0; i < nb; ++i) {
    b[3 * i + 0] = points2[i][0];
    b[3 * i + 1] = points2[i][1];
    b[3 * i + 2] = points2[i][2];
  }

  auto result = py::array_t<T>({na, nb});
  auto buf = result.request();
  T* ptr = static_cast<T*>(buf.ptr);

  distopia::DistanceArrayNoBox(a.data(), b.data(), na, nb, ptr);

  return result;
}

template <typename T>
static py::array_t<T> distance(const std::vector<std::vector<T>> &points) {
  return distance(points, points);
}
} // namespace

namespace lahuta::bindings {
void bind_distance(py::module_ &m) {

  py::class_<lahuta::ContiguousMatrix<double>>(m, "ContiguousMatrix", py::buffer_protocol())
    .def(py::init<int, int>(), py::arg("rows"), py::arg("cols"))
    .def("rows", &lahuta::ContiguousMatrix<double>::rows)
    .def("cols", &lahuta::ContiguousMatrix<double>::cols)

    .def("to_numpy", [](const lahuta::ContiguousMatrix<double>&m) -> py::array_t<double> { return lahuta::numpy::matrix_to_numpy(m); })
    .def("to_numpy_view", [](lahuta::ContiguousMatrix<double> &m) -> py::array_t<double> {
      // zero-copy view tied to the lifetime of this C++ instance
      return py::array_t<double>({m.rows(), m.cols()},
                                 {m.cols() * (py::ssize_t)sizeof(double), (py::ssize_t)sizeof(double)},
                                 m.data(), py::cast(m));
    })
    .def("__call__", [](lahuta::ContiguousMatrix<double>&m, int i, int j) -> double { return m(i, j); }, py::arg("i"), py::arg("j"))
    .def("__getitem__", [](lahuta::ContiguousMatrix<double> &m, py::tuple index) -> double& {
      if (index.size() != 2) 
        throw py::index_error("ContiguousMatrix index must be a 2-tuple");
      int i = py::cast<int>(index[0]);
      int j = py::cast<int>(index[1]);
      return m(i, j);
    }, py::return_value_policy::reference_internal)
    .def("__setitem__", [](lahuta::ContiguousMatrix<double> &m, py::tuple index, double value) {
      if (index.size() != 2) 
        throw py::index_error("ContiguousMatrix index must be a 2-tuple");
      int i = py::cast<int>(index[0]);
      int j = py::cast<int>(index[1]);
      m(i, j) = value;
    })
    .def("__repr__", [](const lahuta::ContiguousMatrix<double> &m) {
      return "ContiguousMatrix(" + std::to_string(m.rows()) + ", " + std::to_string(m.cols()) + ")";
    })
    .def_buffer([](lahuta::ContiguousMatrix<double> &m) -> py::buffer_info {
      return py::buffer_info(
        m.data(),
        sizeof(double),
        py::format_descriptor<double>::format(),
        2,
        { (py::ssize_t)m.rows(), (py::ssize_t)m.cols() },
        { (py::ssize_t)(m.cols() * sizeof(double)), (py::ssize_t)sizeof(double) }
      );
    });

  py::class_<lahuta::DistanceComputation>(m, "DistanceComputation")
    .def_static("distance", [](const std::vector<double> &p1, const std::vector<double> &p2) -> double {
      return lahuta::DistanceComputation::distance(p1, p2);
    }, py::arg("p1"), py::arg("p2"))
    .def_static("distance", [](const std::vector<std::vector<double>> &points1, const std::vector<std::vector<double>> &points2) -> lahuta::ContiguousMatrix<double> {
      return lahuta::DistanceComputation::distance(points1, points2);
    }, py::arg("points1"), py::arg("points2"))
    .def_static("distance", [](const std::vector<std::vector<double>> &points) -> lahuta::ContiguousMatrix<double> {
      return lahuta::DistanceComputation::distance(points);
    }, py::arg("points"))
    .def_static("distance_matrix", [](py::array_t<double, py::array::c_style | py::array::forcecast> points1_np,
                                      py::array_t<double, py::array::c_style | py::array::forcecast> points2_np) -> py::array_t<double> {
      auto points1 = lahuta::numpy::array2vv_double(points1_np);
      auto points2 = lahuta::numpy::array2vv_double(points2_np);
      auto result_matrix = lahuta::DistanceComputation::distance(points1, points2);
      return lahuta::numpy::matrix_to_numpy(result_matrix);
    }, py::arg("points1"), py::arg("points2"))
    .def_static("distance_matrix", [](py::array_t<double, py::array::c_style | py::array::forcecast> points_np) -> py::array_t<double> {
      auto points = lahuta::numpy::array2vv_double(points_np);
      auto result_matrix = lahuta::DistanceComputation::distance(points);
      return lahuta::numpy::matrix_to_numpy(result_matrix);
    }, py::arg("points"))
    .def_static("search", [](const std::vector<std::vector<double>> &points1, const std::vector<std::vector<double>> &points2, double cutoff) -> NSResults {
      return lahuta::DistanceComputation::search(points1, points2, cutoff);
    }, py::arg("points1"), py::arg("points2"), py::arg("cutoff"))
    .def_static("search", [](const std::vector<std::vector<double>> &points, double cutoff) -> NSResults {
      return lahuta::DistanceComputation::search(points, cutoff);
    }, py::arg("points"), py::arg("cutoff"))
    .def_static("search_numpy", [](py::array_t<double, py::array::c_style | py::array::forcecast> points1_np,
                                   py::array_t<double, py::array::c_style | py::array::forcecast> points2_np,
                                   double cutoff) -> NSResults {
      auto points1 = lahuta::numpy::array2vv_double(points1_np);
      auto points2 = lahuta::numpy::array2vv_double(points2_np);
      return lahuta::DistanceComputation::search(points1, points2, cutoff);
    }, py::arg("points1"), py::arg("points2"), py::arg("cutoff"))
    .def_static("search_numpy", [](py::array_t<double, py::array::c_style | py::array::forcecast> points_np,
                                   double cutoff) -> NSResults {
      auto points = lahuta::numpy::array2vv_double(points_np);
      return lahuta::DistanceComputation::search(points, cutoff);
    }, py::arg("points"), py::arg("cutoff"));

}
} // namespace lahuta::bindings
