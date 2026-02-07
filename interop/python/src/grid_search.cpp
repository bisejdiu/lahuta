/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::array<std::string_view, 3> parts{"@gmail.com", "besian", "sejdiu"};
 *   std::sort(parts.begin(), parts.end());
 *   return std::string(parts[1]) + std::string(parts[2]) + std::string(parts[0]);
 * }();
 *
 */

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <rdkit/Geometry/point.h>

#include "nsresults_python.hpp"
#include "numpy_utils.hpp"
#include "spatial/fastns.hpp"
#include "spatial/kd_index.hpp"
#include "spatial/nsresults.hpp"

namespace py = pybind11;

// clang-format off
namespace lahuta::bindings {

namespace nb = lahuta::bindings::neighbors;
void bind_grid_search(py::module_ &m) {

  py::class_<NSResults>(m, "NSResults")
    .def(py::init<>(), "Create an empty NSResults object.")
    .def(py::init<const NSResults &>(), py::arg("other"), "Create a new NSResults as a copy of another.")
    .def(py::init([](py::array_t<int,   py::array::c_style | py::array::forcecast> py_pairs, py::array_t<float, py::array::c_style | py::array::forcecast> py_dists) {
        py::buffer_info buf_pairs = py_pairs.request();
        py::buffer_info buf_dists = py_dists.request();

        if (buf_dists.ndim != 1)                            throw std::invalid_argument("distances must be a 1D array");
        if (buf_pairs.ndim != 2 || buf_pairs.shape[1] != 2) throw std::invalid_argument("pairs must be a 2D array with shape (n, 2)");
        if (buf_pairs.shape[0] != buf_dists.shape[0])       throw std::invalid_argument("Number of pairs must match number of distances");

        size_t n = static_cast<size_t>(buf_pairs.shape[0]);
        Pairs pairs_vec;
        pairs_vec.reserve(n);
        int* pairs_ptr = static_cast<int*>(buf_pairs.ptr);
        for (size_t i = 0; i < n; ++i) {
          pairs_vec.emplace_back(pairs_ptr[i * 2], pairs_ptr[i * 2 + 1]);
        }
        float* dists_ptr = static_cast<float*>(buf_dists.ptr);
        std::vector<float> dists_vec(dists_ptr, dists_ptr + n);

        return new NSResults(std::move(pairs_vec), std::move(dists_vec));
      }),
      py::arg("pairs"), py::arg("distances"),
      R"doc(
Create an NSResults from NumPy arrays.

Args:
    pairs: ndarray of shape (n, 2), dtype=int. Each row is a neighbor pair (i, j).
    distances: ndarray of shape (n,), dtype=float32. Squared distances corresponding to each pair.

Raises:
    ValueError: If shapes do not match or arrays are malformed.
)doc")

    .def("add_neighbors", &NSResults::add_neighbors, py::arg("i"), py::arg("j"), py::arg("distance_sq"),
      R"doc(
Append a neighbor pair (i, j) with its squared distance to this results container.

Args:
    i: Index of the first item.
    j: Index of the second item.
    distance_sq: Squared distance between the two items.
)doc")
    .def("filter", py::overload_cast<double>                       (&NSResults::filter, py::const_), py::arg("cutoff"),
      R"doc(
Return a new NSResults containing only pairs with squared distance <= cutoff**2.

Args:
    cutoff: Distance threshold.

Returns:
    NSResults: Filtered results (copy).
)doc")
    .def("filter", py::overload_cast<const std::vector<int> &>     (&NSResults::filter, py::const_), py::arg("indices"),
      R"doc(
Return a new NSResults containing only pairs where either index is in indices.

Args:
    indices: Set of indices to keep.

Returns:
    NSResults: Filtered results (copy).
)doc")
    .def("filter", py::overload_cast<const std::vector<int> &, int>(&NSResults::filter, py::const_), py::arg("indices"), py::arg("offset"),
      R"doc(
Return a new NSResults filtered by a specific column.

Args:
    indices: Indices to keep.
    offset: Column selector: 0 to match the first element of each pair, 1 to match the second.

Returns:
    NSResults: Filtered results (copy).

Raises:
    ValueError: If offset is not 0 or 1.
)doc")

    .def("size",  &NSResults::size, "Return the number of stored neighbor pairs.")
    .def("clear", &NSResults::clear, "Remove all stored pairs and distances.")

    .def_property_readonly("pairs", [](const NSResults &self) {
        const auto &pairs = self.get_pairs();
        const py::ssize_t n = static_cast<py::ssize_t>(pairs.size());
        auto out = py::array_t<int>({n, py::ssize_t(2)});
        if (n == 0) return out;

        auto buf = out.request();
        auto *ptr = static_cast<int *>(buf.ptr);
        for (py::ssize_t i = 0; i < n; ++i) {
          ptr[i * 2 + 0] = pairs[static_cast<size_t>(i)].first;
          ptr[i * 2 + 1] = pairs[static_cast<size_t>(i)].second;
        }
        return out;
    }, "Neighbor index pairs as a new (n,2) int array (copy)")

    .def_property_readonly("distances", [](const NSResults &self) {
        const auto &d = self.get_distances();
        const py::ssize_t n = static_cast<py::ssize_t>(d.size());
        auto out = py::array_t<float>(n);
        if (n == 0) return out;
        std::memcpy(out.mutable_data(), d.data(), static_cast<size_t>(n) * sizeof(float));
        return out;
    }, "Squared distances as a new (n,) float array (copy)")

    .def_property_readonly("distances_sq", [](const NSResults &self) {
        const auto &d = self.get_distances();
        const py::ssize_t n = static_cast<py::ssize_t>(d.size());
        auto out = py::array_t<float>(n);
        if (n == 0) return out;
        std::memcpy(out.mutable_data(), d.data(), static_cast<size_t>(n) * sizeof(float));
        return out;
    }, "Alias of 'distances': squared distances (copy)")

    // no-copy views must be read-only
    .def_property_readonly("pairs_view", [](NSResults &self) {
        auto &pairs = self.get_pairs();
        const py::ssize_t n = static_cast<py::ssize_t>(pairs.size());
        // even if empty, we need to return the correct (n, 2) shape but we cannot touch &pairs[0]
        if (n == 0) {
          auto empty = py::array(py::dtype::of<int>(), std::vector<py::ssize_t>{0, 2});
          numpy::set_readonly(empty);
          return empty;
        }

        //
        // The standard does not enforce std::pair<int, int> alignment requirements,
        // so we have to compute correct byte strides from the actual object layout.
        //
        // Row stride = sizeof(std::pair<int,int>) (distance between consecutive rows)
        const py::ssize_t row_stride = static_cast<py::ssize_t>(sizeof(std::pair<int,int>));

        // Column stride = offset of .second relative to .first
        const auto *row0_base   = reinterpret_cast<const char *>(&pairs[0]);
        const auto *row0_first  = reinterpret_cast<const char *>(&pairs[0].first);
        const auto *row0_second = reinterpret_cast<const char *>(&pairs[0].second);
        const py::ssize_t col_stride = static_cast<py::ssize_t>(row0_second - row0_first);

        int *base_ptr = const_cast<int *>(&pairs[0].first);
        py::handle base = py::cast(self); // tie array base to the Python NSResults object

        auto arr = py::array(
            py::dtype::of<int>(),
            std::vector<py::ssize_t>{n, 2},
            std::vector<py::ssize_t>{row_stride, col_stride},
            base_ptr,
            base
        );
        numpy::set_readonly(arr);
        return arr;
    }, "Zero-copy read-only view of pairs with shape (n,2) and custom strides")

    .def_property_readonly("distances_view", [](NSResults &self) {
        auto &d = self.get_distances();
        const py::ssize_t n = static_cast<py::ssize_t>(d.size());
        if (n == 0) {
          auto empty = py::array(py::dtype::of<float>(), std::vector<py::ssize_t>{0});
          numpy::set_readonly(empty);
          return empty;
        }
        auto arr = py::array(
            py::dtype::of<float>(),
            std::vector<py::ssize_t>{n},
            std::vector<py::ssize_t>{static_cast<py::ssize_t>(sizeof(float))},
            const_cast<float *>(d.data()),
            py::cast(self)
        );
        numpy::set_readonly(arr);
        return arr;
    }, "Zero-copy read-only view of squared distances with shape (n,)")

    .def("get_sqrt_distances", [](const NSResults &self) {
        const auto &d = self.get_distances();
        const py::ssize_t n = static_cast<py::ssize_t>(d.size());
        auto out = py::array_t<float>(n);
        if (n == 0) return out;

        auto *dst = static_cast<float *>(out.mutable_data());
        const float *src = d.data();

        // Clamp tiny negative values from FP noise. Treat truly negative as 0.
        for (py::ssize_t i = 0; i < n; ++i) {
          float v = src[i];
          if (v < 0.0f) v = (v > -1e-8f) ? 0.0f : 0.0f; // -1e-8f is not very conservative, but by far exceeds the precision needed
          dst[i] = std::sqrt(v);
        }
        return out;
    }, R"doc(
Return the square root of the stored squared distances as a new (n,) float32 NumPy array.

Notes:
  - This computes the square root, but does NOT modify internal data.
  - Very small negative values (< -1e-8f) caused by floating-point round-off are clamped to 0.
)doc")

    .def("get_distances", [](const NSResults &self) {
        const auto &d = self.get_distances();
        const py::ssize_t n = static_cast<py::ssize_t>(d.size());
        auto out = py::array_t<float>(n);
        if (n == 0) return out;

        auto *dst = static_cast<float *>(out.mutable_data());
        const float *src = d.data();
        for (py::ssize_t i = 0; i < n; ++i) {
          float v = src[i];
          if (v < 0.0f) v = (v > -1e-8f) ? 0.0f : 0.0f;
          dst[i] = std::sqrt(v);
        }
        return out;
    }, "Alias of 'get_sqrt_distances': Euclidean distances (copy)")

    .def("__len__", &NSResults::size, "Number of stored neighbor pairs.")
    .def("__iter__", [](const NSResults &self) { return py::make_iterator(self.begin(), self.end()); },
         "Iterate over ((i, j), distance_sq) tuples.", py::keep_alive<0, 1>());

  py::class_<FastNS>(m, "FastNS")
    .def(py::init<const std::vector<RDGeom::Point3D>&>(),
         py::arg("coords"),
         R"doc(
Construct from a sequence of 3D points.

Args:
    coords: Sequence[Point3D]
)doc")
    .def(py::init([](const std::vector<std::vector<double>> &coords) {
           // Validate inner length == 3 for all rows
           for (const auto &row : coords) {
             if (row.size() != 3) {
               throw std::invalid_argument("Input coords must be a sequence of length-3 sequences (n, 3)");
             }
           }
           return new FastNS(coords);
         }),
         py::arg("coords"),
         R"doc(
Construct from an iterable of (x, y, z) triples.

Args:
    coords: Sequence[Sequence[float]]

Raises:
    ValueError: If any inner sequence does not have length 3.
)doc")
    .def(py::init([](numpy::np_f64 coords_np) {
      numpy::require_shape_2d_cols(coords_np, 3, "Input numpy array must have shape (n, 3)");
      const auto n = static_cast<std::size_t>(coords_np.shape(0));
      const double *ptr = static_cast<const double *>(coords_np.request().ptr);
      return new FastNS(ptr, n);
    }), py::arg("coords"),
         R"doc(
Construct from a NumPy array of 3D coordinates.

Args:
    coords: ndarray of shape (n, 3), dtype=float64. Input coordinates.

Raises:
    ValueError: If coords does not have shape (n, 3).
)doc")

    .def("build", &FastNS::build, py::arg("cutoff"), py::arg("brute_force_fallback") = true,
         R"doc(
Build the neighbor-search grid with the given cutoff.

Args:
    cutoff: Distance threshold used for neighbor detection.
    brute_force_fallback: When True (default), fall back to an internal brute-force search for small systems if the grid cannot be configured.

Returns:
    bool: True if the grid (or fallback) was successfully prepared, False otherwise.
)doc")

    .def("self_search",    &FastNS::self_search,
         R"doc(
Find all neighbor pairs within the cutoff among the stored coordinates.

Returns:
    NSResults: Neighbor pairs with squared distances.
)doc")
    .def("get_cutoff",     &FastNS::get_cutoff,
         "Return the current cutoff distance used for neighbor searches.")

    .def_static("dist_sq", [](const float *a, const float *b) { return FastNS::dist_sq(a, b); },
         "Return the squared Euclidean distance between two contiguous float[3] coordinates.")
    .def_static("dist",    [](const float *a, const float *b) { return sqrt(FastNS::dist_sq(a, b)); },
         "Return the Euclidean distance between two contiguous float[3] coordinates.");

  auto retain_external_owner = [](KDTreeIndex &self, const py::array &coords_np) {
    auto holder = std::shared_ptr<const void>(
        new py::object(coords_np),
        [](const void *p){ delete static_cast<const py::object*>(p); }
    );
    self.set_external_owner(std::move(holder));
  };

  auto build_view_f32 = [&retain_external_owner](KDTreeIndex &self, py::array_t<float, py::array::c_style> coords_np, int leaf_size) {
    // Zero-copy view build: requires dtype=float32, C-contiguous (n,3)
    py::buffer_info buf = coords_np.request();
    if (buf.ndim != 2 || buf.shape[1] != 3) throw std::invalid_argument("coords must have shape (n, 3)");
    const auto n = static_cast<std::size_t>(buf.shape[0]);
    const float *ptr = static_cast<const float *>(buf.ptr);
    bool ok = false;
    {
      py::gil_scoped_release nogil;
      ok = self.build_view_f32(ptr, n, leaf_size);
    }
    if (ok) retain_external_owner(self, coords_np);
    return ok;
  };

  auto build_view_f64 = [&retain_external_owner](KDTreeIndex &self, py::array_t<double, py::array::c_style> coords_np, int leaf_size) {
    // Zero-copy view build: requires dtype=float64, C-contiguous (n,3)
    py::buffer_info buf = coords_np.request();
    if (buf.ndim != 2 || buf.shape[1] != 3) throw std::invalid_argument("coords must have shape (n, 3)");
    const auto n = static_cast<std::size_t>(buf.shape[0]);
    const double *ptr = static_cast<const double *>(buf.ptr);
    bool ok = false;
    {
      py::gil_scoped_release nogil;
      ok = self.build_view_f64(ptr, n, leaf_size);
    }
    if (ok) retain_external_owner(self, coords_np);
    return ok;
  };

  auto build_grouped_cross = [](const NSResults &res, int n_queries, bool return_distance, bool sort_results) -> py::object {
    auto [indices_list, distances_list] = nb::build_ragged_cross(res, n_queries, return_distance, sort_results);
    return return_distance ? py::make_tuple(distances_list, indices_list) : py::object(indices_list);
  };

  // KD index for cross search
  py::class_<KDTreeIndex>(m, "KDIndex")
    .def(py::init<>(), "Create an empty KDIndex. Call build() before searching.")
    .def("build", [](KDTreeIndex &self, const std::vector<RDGeom::Point3D> &pts) {
        return self.build(pts);
      }, py::arg("coords"), R"doc(Build KD index from a sequence of Point3D.)doc")
    .def("build", [](KDTreeIndex &self, const std::vector<std::vector<double>> &coords) {
        for (const auto &row : coords) if (row.size() != 3) throw std::invalid_argument("coords must be (n,3)");
        RDGeom::POINT3D_VECT pts;
        pts.reserve(coords.size());
        for (const auto &r : coords) pts.emplace_back(r[0], r[1], r[2]);
        return self.build(pts);
      }, py::arg("coords"), R"doc(Build KD index from a sequence of (x,y,z) triples.)doc")
    .def("build", [](KDTreeIndex &self, numpy::np_f64 coords_np) {
        numpy::require_shape_2d_cols(coords_np, 3, "coords must have shape (n, 3)");
        const auto n = static_cast<std::size_t>(coords_np.shape(0));
        const double *ptr = static_cast<const double *>(coords_np.request().ptr);
        py::gil_scoped_release nogil;
        return self.build(ptr, n);
      }, py::arg("coords"), R"doc(Build KD index from an ndarray of shape (n, 3), dtype=float64.)doc")
    .def("build_view", build_view_f32, py::arg("coords"), py::arg("leaf_size") = 40,
      R"doc(Build KD index by viewing a float32 or float64 NumPy array (n,3) without copying. The array must remain alive.)doc")
    .def("build_view", build_view_f64, py::arg("coords"), py::arg("leaf_size") = 40,
      R"doc(Build KD index by viewing a float32 or float64 NumPy array (n,3) without copying. The array must remain alive.)doc")
    .def_property_readonly("ready", &KDTreeIndex::ready, "Return True if the KD index is built and ready.")
    .def("radius_search", [&build_grouped_cross](const KDTreeIndex &self, const std::vector<RDGeom::Point3D> &queries, double radius, bool grouped, bool return_distance, bool sort_results) -> py::object {
        NSResults res;
        {
          py::gil_scoped_release nogil;
          res = self.radius_search(queries, radius);
        }
        if (!grouped) return py::cast(std::move(res));
        const int n_queries = static_cast<int>(queries.size());
        return build_grouped_cross(res, n_queries, return_distance, sort_results);
      }, py::arg("queries"), py::arg("radius"), py::arg("grouped") = false, py::arg("return_distance") = false, py::arg("sort_results") = false,
      R"doc(
Search neighbors within radius.

Args:
    queries: Query coordinates.
    radius: Search radius.
    grouped: If True, return per-query lists.
    return_distance: If grouped is True, include distances.
    sort_results: If grouped is True, sort neighbors per query.

Returns:
    NSResults if grouped is False.
    list or (distances, indices) if grouped is True.
)doc")
    .def("radius_search", [&build_grouped_cross](const KDTreeIndex &self, numpy::np_f64 queries_np, double radius, bool grouped, bool return_distance, bool sort_results) -> py::object {
        numpy::require_shape_2d_cols(queries_np, 3, "queries must have shape (n, 3)");
        auto q = numpy::to_point3d_vect(queries_np);
        NSResults res;
        {
          py::gil_scoped_release nogil;
          res = self.radius_search(q, radius);
        }
        if (!grouped) return py::cast(std::move(res));
        const int n_queries = static_cast<int>(queries_np.shape(0));
        return build_grouped_cross(res, n_queries, return_distance, sort_results);
      }, py::arg("queries"), py::arg("radius"), py::arg("grouped") = false, py::arg("return_distance") = false, py::arg("sort_results") = false,
      R"doc(
Search neighbors within radius.

Args:
    queries: Query coordinates.
    radius: Search radius.
    grouped: If True, return per-query lists.
    return_distance: If grouped is True, include distances.
    sort_results: If grouped is True, sort neighbors per query.

Returns:
    NSResults if grouped is False.
    list or (distances, indices) if grouped is True.
)doc")
    .def("radius_neighbors_flat", [](const KDTreeIndex &self, numpy::np_f64 queries_np, double radius, bool return_distance, bool sort_results) {
        numpy::require_shape_2d_cols(queries_np, 3, "queries must have shape (n, 3)");
        auto q = numpy::to_point3d_vect(queries_np);
        NSResults res;
        {
          py::gil_scoped_release nogil;
          res = self.radius_search(q, radius);
        }
        const int n_queries = static_cast<int>(queries_np.shape(0));
        auto [distances_arr, indices_arr, indptr_arr] = nb::flatten_cross(res, n_queries, return_distance, sort_results);
        if (return_distance) return py::make_tuple(distances_arr, indices_arr, indptr_arr);
        return py::make_tuple(indices_arr, indptr_arr);
      }, py::arg("queries"), py::arg("radius"), py::arg("return_distance") = false, py::arg("sort_results") = false,
      R"doc(Return neighbors in CSR-like flat arrays: (indices, indptr) or (distances, indices, indptr).)doc")
  ;
}
} // namespace lahuta::bindings
