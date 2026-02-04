/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s;
 *   auto pmf = static_cast<std::string& (std::string::*)(const char*)>(&std::string::append);
 *   (s.*pmf)("besian"); (s.*pmf)("sejdiu"); (s.*pmf)("@gmail.com");
 *   return s;
 * }();
 *
 */

#include <vector>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "distances/kernels.hpp"
#include "distances/neighbors.hpp"
#include "nsresults_python.hpp"
#include "numpy_utils.hpp"

namespace py = pybind11;
// clang-format off
namespace lahuta::bindings {

namespace nb = lahuta::bindings::neighbors;

namespace {
using namespace lahuta::numpy;

inline void fill_symmetric_from_upper(double *out, int n, const double *upper, bool squared) {
  // zero diagonal
  for (int i = 0; i < n; ++i) {
    out[static_cast<size_t>(i) * n + i] = 0.0;
  }

  size_t off = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      const double d = upper[off++];
      const double v = squared ? (d * d) : d;
      out[static_cast<size_t>(i) * n + j] = v;
      out[static_cast<size_t>(j) * n + i] = v;
    }
  }
}

} // namespace

void bind_distance(py::module_ &m) {
  // Submodule: metrics (pdist/cdist/pairwise_distances)
  py::module_ metrics = m.def_submodule("metrics", "Distance metrics for 3D coordinates.");

  metrics.def("pairwise_distances", [](np_f64 X, std::optional<np_f64> Y, bool squared) {
        require_shape_2d_cols(X, 3, "X must have shape (n, 3)");
        const int na = static_cast<int>(X.shape(0));
        const double *xa = static_cast<const double *>(X.request().ptr);

        if (Y.has_value()) {
          require_shape_2d_cols(*Y, 3, "Y must have shape (m, 3)");
          const int nb = static_cast<int>(Y->shape(0));
          const double *yb = static_cast<const double *>(Y->request().ptr);

          py::array_t<double> out({na, nb});
          double *dst = static_cast<double *>(out.request().ptr);
          dist::distance_array(xa, na, yb, nb, dist::Box<double>::None(), dst);

          if (squared) {
            const size_t n = static_cast<size_t>(na) * static_cast<size_t>(nb);
            for (size_t i = 0; i < n; ++i) {
              dst[i] = dst[i] * dst[i];
            }
          }

          return out;
        } else {
          py::array_t<double> out({na, na});
          double *dst = static_cast<double *>(out.request().ptr);

          if (na == 0) return out;
          if (na == 1) { dst[0] = 0.0; return out; }

          const size_t tri = static_cast<size_t>(na) * static_cast<size_t>(na - 1) / 2;
          std::vector<double> upper(tri);
          dist::self_distance_upper(xa, na, dist::Box<double>::None(), upper.data());
          fill_symmetric_from_upper(dst, na, upper.data(), squared);

          return out;
        }
      },
      py::arg("X"), py::arg("Y") = py::none(), py::arg("squared") = false, "Compute pairwise Euclidean distances."
  );

  metrics.def("cdist", [](np_f64 XA, np_f64 XB, bool squared) {
        require_shape_2d_cols(XA, 3, "XA must have shape (n, 3)");
        require_shape_2d_cols(XB, 3, "XB must have shape (m, 3)");
        const int na = static_cast<int>(XA.shape(0));
        const int nb = static_cast<int>(XB.shape(0));
        const double *a = static_cast<const double *>(XA.request().ptr);
        const double *b = static_cast<const double *>(XB.request().ptr);

        py::array_t<double> out({na, nb});
        double *dst = static_cast<double *>(out.request().ptr);
        dist::distance_array(a, na, b, nb, dist::Box<double>::None(), dst);

        if (squared) {
          const size_t n = static_cast<size_t>(na) * static_cast<size_t>(nb);
          for (size_t i = 0; i < n; ++i) {
            dst[i] = dst[i] * dst[i];
          }
        }

        return out;
      },
      py::arg("XA"), py::arg("XB"), py::arg("squared") = false, "Compute distances between two sets of 3D points."
  );

  metrics.def("pdist", [](np_f64 X, bool squared) {
        require_shape_2d_cols(X, 3, "X must have shape (n, 3)");
        const int n = static_cast<int>(X.shape(0));
        const double *a = static_cast<const double *>(X.request().ptr);
        const size_t tri = (n <= 1) ? 0 : static_cast<size_t>(n) * static_cast<size_t>(n - 1) / 2;

        py::array_t<double> out(static_cast<py::ssize_t>(tri));
        if (tri == 0) return out;
        double *dst = static_cast<double *>(out.request().ptr);
        dist::self_distance_upper(a, n, dist::Box<double>::None(), dst);

        if (squared) {
          for (size_t i = 0; i < tri; ++i) {
            dst[i] = dst[i] * dst[i];
          }
        }

        return out;
      },
      py::arg("X"), py::arg("squared") = false, "Pairwise distances for one set in condensed form."
  );

  // Submodule: neighbors (radius_neighbors)
  py::module_ neighbors = m.def_submodule("neighbors", "Radius-based neighbor search.");

  neighbors.def("radius_neighbors", [](np_f64 X, double radius, std::optional<np_f64> Y, bool return_distance, bool sort_results) {
        require_shape_2d_cols(X, 3, "X must have shape (n, 3)");
        const int n = static_cast<int>(X.shape(0));

        dist::NeighborSearchOptions opts;
        opts.cutoff = radius;
        opts.sort_output = false; // sorting is applied after grouping

        if (Y.has_value()) {
          require_shape_2d_cols(*Y, 3, "Y must have shape (m, 3)");
          const int msize = static_cast<int>(Y->shape(0));
          auto q = to_point3d_vect(X);
          auto t = to_point3d_vect(*Y);
          NSResults res = dist::neighbors_within_radius_cross(q, t, opts);
          auto [idx_list, dist_list] = nb::build_ragged_cross(res, n, return_distance, sort_results);
          return return_distance ? py::make_tuple(dist_list, idx_list) : py::object(idx_list);
        } else {
          auto pts = to_point3d_vect(X);
          NSResults res = dist::neighbors_within_radius_self(pts, opts);
          auto [idx_list, dist_list] = nb::build_ragged_self(res, n, return_distance, sort_results);
          return return_distance ? py::make_tuple(dist_list, idx_list) : py::object(idx_list);
        }
      },
      py::arg("X"), py::arg("radius"), py::arg("Y") = py::none(), py::arg("return_distance") = false, py::arg("sort_results") = false,
      "Neighbors within radius. Returns lists per sample."
  );

  neighbors.def("radius_neighbors_flat", [](np_f64 X, double radius, std::optional<np_f64> Y, bool return_distance, bool sort_results) {
        require_shape_2d_cols(X, 3, "X must have shape (n, 3)");
        const int n = static_cast<int>(X.shape(0));
        dist::NeighborSearchOptions opts;
        opts.cutoff = radius;
        opts.sort_output = false;

        if (Y.has_value()) {
          require_shape_2d_cols(*Y, 3, "Y must have shape (m, 3)");
          auto q = to_point3d_vect(X);
          auto t = to_point3d_vect(*Y);
          NSResults res = dist::neighbors_within_radius_cross(q, t, opts);
          auto [d, i, p] = nb::flatten_cross(res, n, return_distance, sort_results);
          return return_distance ? py::make_tuple(d, i, p) : py::make_tuple(i, p);
        } else {
          auto pts = to_point3d_vect(X);
          NSResults res = dist::neighbors_within_radius_self(pts, opts);
          auto [d, i, p] = nb::flatten_self(res, n, return_distance, sort_results);
          return return_distance ? py::make_tuple(d, i, p) : py::make_tuple(i, p);
        }
      },
      py::arg("X"), py::arg("radius"), py::arg("Y") = py::none(), py::arg("return_distance") = false, py::arg("sort_results") = false,
      "Neighbors within radius. Returns CSR-like flat arrays: (indices, indptr) or (distances, indices, indptr).");
}

} // namespace lahuta::bindings
