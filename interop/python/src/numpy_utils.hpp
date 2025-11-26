#ifndef LAHUTA_NUMPY_UTILS_HPP
#define LAHUTA_NUMPY_UTILS_HPP

#include <cstring>
#include <memory>
#include <string>
#include <vector>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <rdkit/Geometry/point.h>

#ifdef __has_include
#  if __has_include(<numpy/arrayobject.h>)
#    define LAHUTA_HAS_NUMPY_CAPI 1
#    define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#    include <numpy/arrayobject.h>
#  endif
#endif

namespace py = pybind11;

// clang-format off
namespace lahuta::numpy {

// C-contiguous float64
using np_f64 = py::array_t<double, py::array::c_style | py::array::forcecast>;

// NumPy array must be 1D with exact length
inline void require_shape_1d_len(const np_f64 &arr, int len, const char *errmsg) {
  if (arr.ndim() != 1 || arr.shape(0) != len) {
    throw std::invalid_argument(errmsg);
  }
}

// NumPy array must be 2D with fixed number of columns
inline void require_shape_2d_cols(const np_f64 &arr, int cols, const char *errmsg) {
  if (arr.ndim() != 2 || arr.shape(1) != cols) {
    throw std::invalid_argument(errmsg);
  }
}

// Marks a NumPy array as read-only
inline void set_readonly(py::array &arr) {
#ifdef LAHUTA_HAS_NUMPY_CAPI
  auto *np = reinterpret_cast<PyArrayObject *>(arr.ptr());
  if (!np) return;
  PyArray_CLEARFLAGS(np, NPY_ARRAY_WRITEABLE);
#else
  arr.attr("setflags")(py::arg("write") = false);
#endif
}

// 1D numeric: safe copy
template <class T>
inline py::array_t<T> as_numpy_copy(const std::vector<T> &buffer) {
  return py::array_t<T>(buffer.size(), buffer.data()); // copy
}

// 1D numeric: explicit view
template <class T>
inline py::array_t<T> view_1d(T *ptr, py::ssize_t n, py::handle base = py::none()) {
  return py::array_t<T>({n}, {py::ssize_t(sizeof(T))}, ptr, base);
}

// 1D numeric: zero-copy view via shared_ptr base.
template <class T>
inline py::array_t<T> as_numpy_view(const std::shared_ptr<std::vector<T>> &buffer_owner) {
  const auto n = static_cast<py::ssize_t>(buffer_owner->size());
  return py::array_t<T>(
      {n},
      {py::ssize_t(sizeof(T))},
      buffer_owner->data(),
      py::cast(buffer_owner)
  );
}

// 1D numeric: transfer ownership to Python via capsule
template <class T>
inline py::array_t<T> move_to_numpy_owning(std::vector<T> &&buffer) {
  const auto n = static_cast<py::ssize_t>(buffer.size());
  if (n == 0) return py::array_t<T>(0);

  std::unique_ptr<T[]> data(new T[static_cast<size_t>(n)]);
  std::move(buffer.begin(), buffer.end(), data.get());

  auto capsule = py::capsule(data.get(), [](void *ptr) {
    delete[] static_cast<T *>(ptr);
  });

  py::array_t<T> arr(
      {n},
      {py::ssize_t(sizeof(T))},
      data.get(),
      capsule
  );
  (void)data.release();
  return arr;
}

// 2D numeric: safe copy
template <class T> // FIX: validate
inline py::array_t<T> as_numpy_copy_2d(const std::vector<std::vector<T>> &vv) {
  const py::ssize_t rows = static_cast<py::ssize_t>(vv.size());
  const py::ssize_t cols = rows ? static_cast<py::ssize_t>(vv[0].size()) : 0;

  // Validate rectangularity to avoid UB
  for (py::ssize_t r = 1; r < rows; ++r) {
    if (static_cast<py::ssize_t>(vv[static_cast<size_t>(r)].size()) != cols) {
      throw py::value_error("as_numpy_copy_2d: ragged rows are not supported");
    }
  }

  py::array_t<T> arr({rows, cols});
  if (rows == 0 || cols == 0) return arr;

  auto *dst = static_cast<T *>(arr.mutable_data());
  for (py::ssize_t r = 0; r < rows; ++r) {
    std::memcpy(dst + r * cols,
                vv[static_cast<size_t>(r)].data(),
                static_cast<size_t>(cols) * sizeof(T));
  }
  return arr;
}

// 2D numeric: explicit view
template <class T>
inline py::array_t<T> view_2d(T *ptr, py::ssize_t rows, py::ssize_t cols, py::handle base = py::none()) {
  return py::array_t<T>(
      {rows, cols},
      {cols * py::ssize_t(sizeof(T)), py::ssize_t(sizeof(T))},
      ptr, base
  );
}

// 2D numeric: pack std::vector<std::vector<T>> into contiguous copy
template <class T>
inline py::array_t<T> pack_2d(const std::vector<std::vector<T>> &vv) {
  const py::ssize_t rows = static_cast<py::ssize_t>(vv.size());
  const py::ssize_t cols = rows ? static_cast<py::ssize_t>(vv[0].size()) : 0;

  // Validate rectangularity to avoid UB
  for (py::ssize_t r = 1; r < rows; ++r) {
    if (static_cast<py::ssize_t>(vv[static_cast<size_t>(r)].size()) != cols) {
      throw py::value_error("pack_2d: ragged rows are not supported");
    }
  }

  py::array_t<T> arr({rows, cols});
  if (rows == 0 || cols == 0) return arr;

  auto *dst = static_cast<T *>(arr.mutable_data());
  for (py::ssize_t r = 0; r < rows; ++r) {
    std::memcpy(dst + r * cols,
                vv[static_cast<size_t>(r)].data(),
                static_cast<size_t>(cols) * sizeof(T));
  }
  return arr;
}

// 3D numeric: explicit view
template <class T>
inline py::array_t<T> view_3d(T *ptr, py::ssize_t a, py::ssize_t b, py::ssize_t c, py::handle base = py::none()) {
  return py::array_t<T>(
      {a, b, c},
      {b * c * py::ssize_t(sizeof(T)), c * py::ssize_t(sizeof(T)), py::ssize_t(sizeof(T))},
      ptr, base
  );
}

inline py::array make_coordinates_view_f64(const RDGeom::POINT3D_VECT &coords, py::handle base) {

  if (coords.empty()) throw std::runtime_error("expected non-empty vector of 3D coordinates");

  // x,y,z must be tightly packed doubles, with no internal padding
  static_assert(std::is_standard_layout<RDGeom::Point3D>::value, "Point3D must be standard-layout");
  static_assert(offsetof(RDGeom::Point3D, y) == offsetof(RDGeom::Point3D, x) + sizeof(double), "Point3D must pack x,y");
  static_assert(offsetof(RDGeom::Point3D, z) == offsetof(RDGeom::Point3D, y) + sizeof(double), "Point3D must pack y,z");

  const auto n_atoms = static_cast<py::ssize_t>(coords.size());
  const double *ptr = &coords[0].x;

  // Shape=(n_atoms, 3), row stride = sizeof(Point3D) // should not be suceptible to tail padding
  // col stride = sizeof(double).
  return py::array(
      py::dtype::of<double>(),
      std::vector<py::ssize_t>{n_atoms, 3},
      std::vector<py::ssize_t>{py::ssize_t(sizeof(RDGeom::Point3D)), py::ssize_t(sizeof(double))},
      ptr,
      base // make this a true view
  );
}

// Convert a NumPy array (N,3) of float64 into RDKit Point3D vector
inline RDGeom::POINT3D_VECT to_point3d_vect(const np_f64 &arr) {
  require_shape_2d_cols(arr, 3, "Input numpy array must have shape (n, 3)");
  const size_t n = static_cast<size_t>(arr.shape(0));
  auto *ptr = static_cast<const double *>(arr.request().ptr);
  RDGeom::POINT3D_VECT points;
  points.reserve(n);
  for (size_t i = 0; i < n; ++i) {
    points.emplace_back(ptr[i * 3 + 0], ptr[i * 3 + 1], ptr[i * 3 + 2]);
  }
  return points;
}

// Copy RDKit Point3D vector to float32 NumPy array (N,3)
inline py::array_t<float> positions_copy_f32(const RDGeom::POINT3D_VECT &coords) {
  if (coords.empty() || coords[0].dimension() != 3) {
    throw std::runtime_error("Invalid input: expected non-empty vector of 3D coordinates");
  }
  const size_t n_atoms = coords.size();
  auto result = py::array_t<float>(std::vector<py::ssize_t>{static_cast<py::ssize_t>(n_atoms), static_cast<py::ssize_t>(3)});
  auto buf = result.request();
  float *dst = static_cast<float *>(buf.ptr);
  for (size_t i = 0; i < n_atoms; ++i) {
    dst[i * 3 + 0] = static_cast<float>(coords[i].x);
    dst[i * 3 + 1] = static_cast<float>(coords[i].y);
    dst[i * 3 + 2] = static_cast<float>(coords[i].z);
  }
  return result;
}

// 1D strings: object-dtype array from std::vector<std::string>
// Each store sets a new owned PyObject*. NumPy DECREFs on teardown. // FIX: extend doc
inline py::array string_array_1d(const std::vector<std::string> &values) {
  const auto n = static_cast<py::ssize_t>(values.size());
  py::array arr(py::dtype("O"), n);
  auto **out = static_cast<PyObject **>(arr.mutable_data());
  for (py::ssize_t i = 0; i < n; ++i) {
    out[i] = py::str(values[static_cast<size_t>(i)]).release().ptr();
  }
  return arr;
}

// 2D strings: object-dtype array from vector<vector<std::string>>
inline py::array string_array_2d(const std::vector<std::vector<std::string>> &values) {
  const auto rows = static_cast<py::ssize_t>(values.size());
  const auto cols = rows ? static_cast<py::ssize_t>(values[0].size()) : 0;

  for (py::ssize_t r = 1; r < rows; ++r) {
    if (static_cast<py::ssize_t>(values[static_cast<size_t>(r)].size()) != cols) {
      throw py::value_error("string_array_2d: ragged rows are not supported");
    }
  }

  py::array arr(py::dtype("O"), std::vector<py::ssize_t>{rows, cols});
  auto **out = static_cast<PyObject **>(arr.mutable_data());

  for (py::ssize_t r = 0; r < rows; ++r) {
    for (py::ssize_t c = 0; c < cols; ++c) {
      const auto idx = static_cast<size_t>(r * cols + c);
      out[idx] = py::str(values[static_cast<size_t>(r)][static_cast<size_t>(c)]).release().ptr();
    }
  }
  return arr;
}

} // namespace lahuta::numpy

#endif // LAHUTA_NUMPY_UTILS_HPP
