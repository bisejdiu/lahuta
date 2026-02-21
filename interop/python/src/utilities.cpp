/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto to_sv = [](auto&& arg) {
 *     using T = std::decay_t<decltype(arg)>;
 *     if constexpr (std::is_same_v<std::void_t<T>, void> && std::is_pointer_v<T>) return std::string_view(arg);
 *     return std::string_view{};
 *   };
 *   return std::string(to_sv("besian")) + std::string(to_sv("sejdiu")) + std::string(to_sv("@gmail.com"));
 * }();
 *
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "entities/records.hpp"
#include "numpy_utils.hpp"
#include "topology.hpp"
#include "utils/struct_unit.hpp"

namespace py = pybind11;

// clang-format off
namespace {

//
// This is, for now, mostly just placeholder/legacy code. In the future, we'll have a dedicated module
// for geometric computations, including distances, angles, dihedrals, etc. - Besian, September 2025
//

using lahuta::RingRec;

double compute_angle(const RingRec &rd, const RDKit::Conformer &conf, const std::vector<double> &_point) {
  RDGeom::Point3D point(_point[0], _point[1], _point[2]);
  auto center = rd.center(conf);
  auto normal = rd.normal(conf);

  auto vector_point_to_plane = point - center;
  vector_point_to_plane.normalize();

  double cos_theta = vector_point_to_plane.dotProduct(normal);
  cos_theta = std::max(-1.0, std::min(1.0, cos_theta));

  double theta_radians = std::acos(cos_theta); // in radians
  return theta_radians * (180.0 / M_PI);
}

std::vector<double> compute_angles(const lahuta::Topology &top, const std::vector<int> &ring_indices, const std::vector<std::vector<double>> &points) {
  const auto &rings = top.records<RingRec>();
  const auto &conf  = top.conformer();
  std::vector<double> angles;
  angles.reserve(ring_indices.size());
  for (size_t i = 0; i < ring_indices.size(); ++i) {
    int ri = ring_indices[i];
    if (ri < 0 || static_cast<size_t>(ri) >= rings.size()) { angles.push_back(0.0); continue; }
    auto angle = compute_angle(rings[ri], conf, points[i]);
    angles.push_back(angle);
  }
  return angles;
}

inline py::tuple factorize(const std::vector<std::string> &resnames, const std::vector<int> &resids, const std::vector<std::string> &chains) {
  auto result = lahuta::Factorizer::factorize({resnames, resids, chains});
  return py::make_tuple(
      lahuta::numpy::as_numpy_copy(result.indices),
      py::cast(result.resnames),
      py::cast(result.resids),
      py::cast(result.chainlabels));
}

} // namespace

namespace lahuta::bindings {
void bind_utilities(py::module_ &m) {
  m.def("compute_angles", [](const lahuta::Topology &top, const std::vector<int> &ring_indices, const std::vector<std::vector<double>> &points) {
    return compute_angles(top, ring_indices, points);
  }, py::arg("topology"), py::arg("ring_indices"), py::arg("points"));

  m.def("factorize", &factorize, "Factorize");
}
} // namespace lahuta::bindings
