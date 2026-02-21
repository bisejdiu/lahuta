/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::ostringstream os; os << "besian" << "sejdiu" << "@gmail.com";
 *   return os.str();
 * }();
 *
 */

#ifndef LAHUTA_ANALYSIS_DSSP_ASA_HPP
#define LAHUTA_ANALYSIS_DSSP_ASA_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include <Geometry/point.h>

#include "analysis/dssp/precompute.hpp"
#include "distances/neighbors.hpp"
#include "utils/math_constants.hpp"

namespace lahuta::analysis {

inline constexpr double DsspRadiusN        = 1.65;
inline constexpr double DsspRadiusCA       = 1.87;
inline constexpr double DsspRadiusC        = 1.76;
inline constexpr double DsspRadiusO        = 1.4;
inline constexpr double DsspRadiusSideAtom = 1.8;
inline constexpr double DsspRadiusWater    = 1.4;

struct DsspBox {
  RDGeom::Point3D min;
  RDGeom::Point3D max;
  RDGeom::Point3D center;
  double radius = 0.0;
};

namespace detail {

inline RDGeom::Point3D make_point(double x, double y, double z) noexcept {
  RDGeom::Point3D p;
  p.x = x;
  p.y = y;
  p.z = z;
  return p;
}

inline void extend_box(DsspBox &box, const RDGeom::Point3D &atom, double radius) noexcept {
  box.min.x = std::min(box.min.x, atom.x - radius);
  box.min.y = std::min(box.min.y, atom.y - radius);
  box.min.z = std::min(box.min.z, atom.z - radius);
  box.max.x = std::max(box.max.x, atom.x + radius);
  box.max.y = std::max(box.max.y, atom.y + radius);
  box.max.z = std::max(box.max.z, atom.z + radius);
}

[[nodiscard]] inline DsspBox build_box(const DsspResidue &res) noexcept {
  const double inf = std::numeric_limits<double>::infinity();
  DsspBox box;
  box.min = make_point(inf, inf, inf);
  box.max = make_point(-inf, -inf, -inf);

  const double expand_n  = DsspRadiusN + 2.0 * DsspRadiusWater;
  const double expand_ca = DsspRadiusCA + 2.0 * DsspRadiusWater;
  const double expand_c  = DsspRadiusC + 2.0 * DsspRadiusWater;
  const double expand_o  = DsspRadiusO + 2.0 * DsspRadiusWater;
  const double expand_sc = DsspRadiusSideAtom + 2.0 * DsspRadiusWater;

  extend_box(box, res.n, expand_n);
  extend_box(box, res.ca, expand_ca);
  extend_box(box, res.c, expand_c);
  extend_box(box, res.o, expand_o);

  for (const auto &sc : res.side_chain) {
    extend_box(box, sc.pos, expand_sc);
  }
  if (res.has_oxt) {
    extend_box(box, res.oxt, expand_sc);
  }

  const double dx = box.max.x - box.min.x;
  const double dy = box.max.y - box.min.y;
  const double dz = box.max.z - box.min.z;

  box.radius = std::max(dx, std::max(dy, dz));
  box.center = make_point((box.min.x + box.max.x) * 0.5,
                          (box.min.y + box.max.y) * 0.5,
                          (box.min.z + box.max.z) * 0.5);
  return box;
}

[[nodiscard]] inline bool atom_intersects_box(const RDGeom::Point3D &atom, double radius,
                                              const DsspBox &box) noexcept {
  return atom.x + radius >= box.min.x && atom.x - radius <= box.max.x && atom.y + radius >= box.min.y &&
         atom.y - radius <= box.max.y && atom.z + radius >= box.min.z && atom.z - radius <= box.max.z;
}

struct Accumulator {
  struct Candidate {
    RDGeom::Point3D location;
    double radius   = 0.0;
    double distance = 0.0;
  };

  std::vector<Candidate> items;

  void add(const RDGeom::Point3D &a, const RDGeom::Point3D &b, double d, double r) noexcept {
    const double dist  = distance_sq(a, b);
    d                 += DsspRadiusWater;
    r                 += DsspRadiusWater;
    double test        = d + r;
    test              *= test;
    if (dist < test && dist > 0.0001) {
      items.emplace_back(Candidate{b - a, r * r, dist});
    }
  }

  void sort() noexcept {
    std::sort(items.begin(), items.end(), [](const Candidate &lhs, const Candidate &rhs) noexcept {
      return lhs.distance < rhs.distance;
    });
  }
};

struct DsspSurfaceDots {
  std::vector<RDGeom::Point3D> points;
  double weight = 0.0;

  [[nodiscard]] static const DsspSurfaceDots &instance() noexcept {
    static const DsspSurfaceDots dots(200);
    return dots;
  }

  explicit DsspSurfaceDots(int n) {
    const double p            = 2.0 * static_cast<double>(n) + 1.0;
    const double golden_ratio = (1.0 + std::sqrt(5.0)) / 2.0;
    weight                    = (4.0 * lahuta::Pi) / p;
    points.reserve(static_cast<std::size_t>(p));
    for (int i = -n; i <= n; ++i) {
      const double lat = std::asin((2.0 * static_cast<double>(i)) / p);
      const double lon = std::fmod(static_cast<double>(i), golden_ratio) * 2.0 * lahuta::Pi / golden_ratio;
      const double cos_lat = std::cos(lat);
      points.emplace_back(std::sin(lon) * cos_lat, std::cos(lon) * cos_lat, std::sin(lat));
    }
  }
};

[[nodiscard]] inline std::vector<std::vector<std::size_t>>
build_neighbor_lists(const std::vector<DsspBox> &boxes) {
  if (boxes.empty()) return {};

  double max_radius = 0.0;
  for (const auto &box : boxes) {
    max_radius = std::max(max_radius, box.radius);
  }
  const double cutoff = 2.0 * max_radius;

  RDGeom::POINT3D_VECT centers;
  centers.reserve(boxes.size());
  for (const auto &box : boxes) {
    centers.emplace_back(box.center);
  }

  dist::NeighborSearchOptions opts;
  opts.cutoff               = cutoff;
  opts.brute_force_fallback = true;
  opts.sort_output          = false;

  const auto results = dist::neighbors_within_radius_self(centers, opts);

  std::vector<std::vector<std::size_t>> neighbors(boxes.size());
  for (auto &n : neighbors) {
    n.reserve(32);
  }

  for (const auto &pair : results.get_pairs()) {
    if (pair.first < 0 || pair.second < 0) continue;
    const auto i = static_cast<std::size_t>(pair.first);
    const auto j = static_cast<std::size_t>(pair.second);
    if (i >= boxes.size() || j >= boxes.size()) continue;

    const double dist2      = distance_sq(boxes[i].center, boxes[j].center);
    const double box_cutoff = boxes[i].radius + boxes[j].radius;
    if (dist2 < box_cutoff * box_cutoff) {
      neighbors[i].push_back(j);
      neighbors[j].push_back(i);
    }
  }

  for (std::size_t i = 0; i < boxes.size(); ++i) {
    neighbors[i].push_back(i);
  }

  return neighbors;
}

[[nodiscard]] inline double calculate_surface_atom(const RDGeom::Point3D &atom, double radius,
                                                   const std::vector<DsspResidue> &residues,
                                                   const std::vector<std::vector<std::size_t>> &neighbors,
                                                   const std::vector<DsspBox> &boxes,
                                                   std::size_t idx) noexcept {
  Accumulator acc;
  for (const auto n_idx : neighbors[idx]) {
    if (n_idx >= residues.size()) continue;
    if (!atom_intersects_box(atom, radius, boxes[n_idx])) continue;

    const auto &res = residues[n_idx];
    acc.add(atom, res.n, radius, DsspRadiusN);
    acc.add(atom, res.ca, radius, DsspRadiusCA);
    acc.add(atom, res.c, radius, DsspRadiusC);
    acc.add(atom, res.o, radius, DsspRadiusO);
    for (const auto &sc : res.side_chain) {
      acc.add(atom, sc.pos, radius, DsspRadiusSideAtom);
    }
    if (res.has_oxt) {
      acc.add(atom, res.oxt, radius, DsspRadiusSideAtom);
    }
  }

  acc.sort();

  const double expanded = radius + DsspRadiusWater;
  double surface        = 0.0;

  const auto &dots = DsspSurfaceDots::instance();
  for (const auto &p : dots.points) {
    const RDGeom::Point3D xx = make_point(p.x * expanded, p.y * expanded, p.z * expanded);
    bool free                = true;
    for (const auto &c : acc.items) {
      if (c.radius >= distance_sq(xx, c.location)) {
        free = false;
        break;
      }
    }
    if (free) surface += dots.weight;
  }

  return surface * expanded * expanded;
}

} // namespace detail

[[nodiscard]] inline std::vector<double> compute_accessibility(const std::vector<DsspResidue> &residues) {
  std::vector<double> out;
  out.assign(residues.size(), 0.0);
  if (residues.empty()) return out;

  std::vector<DsspBox> boxes;
  boxes.reserve(residues.size());
  for (const auto &res : residues) {
    boxes.push_back(detail::build_box(res));
  }

  auto neighbors = detail::build_neighbor_lists(boxes);

  for (std::size_t i = 0; i < residues.size(); ++i) {
    const auto &res = residues[i];
    double surface  = 0.0;

    surface += detail::calculate_surface_atom(res.n, DsspRadiusN, residues, neighbors, boxes, i);
    surface += detail::calculate_surface_atom(res.ca, DsspRadiusCA, residues, neighbors, boxes, i);
    surface += detail::calculate_surface_atom(res.c, DsspRadiusC, residues, neighbors, boxes, i);
    surface += detail::calculate_surface_atom(res.o, DsspRadiusO, residues, neighbors, boxes, i);

    for (const auto &sc : res.side_chain) {
      surface += detail::calculate_surface_atom(sc.pos, DsspRadiusSideAtom, residues, neighbors, boxes, i);
    }
    if (res.has_oxt) {
      surface += detail::calculate_surface_atom(res.oxt, DsspRadiusSideAtom, residues, neighbors, boxes, i);
    }
    out[i] = surface;
  }
  return out;
}

} // namespace lahuta::analysis

#endif // LAHUTA_ANALYSIS_DSSP_ASA_HPP
