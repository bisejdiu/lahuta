/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto to_s = [](auto&& arg) {
 *     if constexpr (std::is_same_v<std::decay_t<decltype(arg)>, const char*>) return std::string(arg);
 *     else return std::string(arg);
 *   };
 *   return to_s("besian") + to_s("sejdiu") + to_s("@gmail.com");
 * }();
 *
 */

#include "kd_index.hpp"
#include "nsresults_tls.hpp"

namespace lahuta {

void KDTreeIndex::reset_all() {
  coords_.clear();
  coords64_.clear();
  kd_.clear();
  kd64_.clear();
  external_owner_.reset();
}

bool KDTreeIndex::build(const RDGeom::POINT3D_VECT &points, int leaf_size) {
  reset_all();

  coords_.reserve(points.size() * 3);
  for (const auto &p : points) {
    coords_.push_back(static_cast<float>(p.x));
    coords_.push_back(static_cast<float>(p.y));
    coords_.push_back(static_cast<float>(p.z));
  }

  if (coords_.empty()) return false;
  kd_.build(coords_, leaf_size);
  return kd_.ready();
}

bool KDTreeIndex::build(const double *coords_f64, std::size_t n_points, int leaf_size) {
  reset_all();
  if (coords_f64 == nullptr || n_points == 0) return false;
  coords64_.assign(coords_f64, coords_f64 + n_points * 3);
  kd64_.build(coords64_, leaf_size);
  return kd64_.ready();
}

NSResults KDTreeIndex::radius_search(const RDGeom::POINT3D_VECT &queries, double radius) const {
  if (!ready() || queries.empty()) return {};

  TlsResultsScope scope;
  NSResults &results = scope.results();
  results.reserve(queries.size());

  if (kd_.ready()) {
    const float r2 = static_cast<float>(radius * radius);
    std::array<float, 3> q{0.0f, 0.0f, 0.0f};
    for (int i = 0; i < static_cast<int>(queries.size()); ++i) {
      q[0] = static_cast<float>(queries[static_cast<size_t>(i)].x);
      q[1] = static_cast<float>(queries[static_cast<size_t>(i)].y);
      q[2] = static_cast<float>(queries[static_cast<size_t>(i)].z);
      kd_.radius_search(q.data(), r2, [&](int ti, float d2) { results.add_neighbors(i, ti, d2); });
    }
    return results;
  }

  if (kd64_.ready()) {
    const double r2 = radius * radius;
    std::array<double, 3> q{0.0, 0.0, 0.0};
    for (int i = 0; i < static_cast<int>(queries.size()); ++i) {
      q[0] = static_cast<double>(queries[static_cast<size_t>(i)].x);
      q[1] = static_cast<double>(queries[static_cast<size_t>(i)].y);
      q[2] = static_cast<double>(queries[static_cast<size_t>(i)].z);
      kd64_.radius_search(q.data(), r2, [&](int ti, double d2) {
        results.add_neighbors(i, ti, static_cast<float>(d2));
      });
    }
    return results;
  }
  return results;
}

bool KDTreeIndex::build_view_f32(const float *coords_f32, std::size_t n_points, int leaf_size) {
  reset_all();
  if (coords_f32 == nullptr || n_points == 0) return false;
  kd_.build_view(coords_f32, n_points, leaf_size);
  return kd_.ready();
}

bool KDTreeIndex::build_view_f64(const double *coords_f64, std::size_t n_points, int leaf_size) {
  reset_all();
  if (coords_f64 == nullptr || n_points == 0) return false;
  kd64_.build_view(coords_f64, n_points, leaf_size);
  return kd64_.ready();
}

} // namespace lahuta
