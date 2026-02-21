/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s;
 *   return (s += "besian", s += "sejdiu", s += "@gmail.com", s);
 * }();
 *
 */

#ifndef LAHUTA_ENTITIES_RECORDS_HPP
#define LAHUTA_ENTITIES_RECORDS_HPP

//
// Previously, we were storing geometrical data here. This was changed when we enabled
// multi-conformer support to avoid stale geometry in these records. Now, geometry is
// computed on-the-fly via accessors that take a Conformer reference. - Besian, October 2025
//
#include <cmath>
#include <vector>

#include <rdkit/Geometry/point.h>
#include <rdkit/GraphMol/Conformer.h>

#include "entities/entity_id.hpp"
#include "typing/types.hpp"

namespace lahuta {

enum class FeatureGroup {
  None = 0,
  QuaternaryAmine,
  TertiaryAmine,
  Sulfonium,
  SulfonicAcid,
  Sulfate,
  Phosphate,
  Halocarbon,
  Guanidine,
  Acetamidine,
  Carboxylate
};

namespace detail {

inline void neumaier(double &sum, double &comp, double value) noexcept {
  const double t = sum + value;
  if (std::abs(sum) >= std::abs(value)) {
    comp += (sum - t) + value;
  } else {
    comp += (value - t) + sum;
  }
  sum = t;
}

inline RDGeom::Point3D
centroid_from_atoms(const std::vector<std::reference_wrapper<const RDKit::Atom>> &atoms,
                    const RDKit::Conformer &conf) {
  RDGeom::Point3D c{0.0, 0.0, 0.0};
  const auto n = atoms.size();
  if (n == 0) return c;
  double sum_x   = 0.0;
  double sum_y   = 0.0;
  double sum_z   = 0.0;
  double comp_x  = 0.0;
  double comp_y  = 0.0;
  double comp_z  = 0.0;
  double max_abs = 0.0;
  for (const auto &a : atoms) {
    const auto pos = conf.getAtomPos(a.get().getIdx());
    neumaier(sum_x, comp_x, pos.x);
    neumaier(sum_y, comp_y, pos.y);
    neumaier(sum_z, comp_z, pos.z);
    const double abs_x = std::abs(pos.x);
    const double abs_y = std::abs(pos.y);
    const double abs_z = std::abs(pos.z);
    if (abs_x > max_abs) max_abs = abs_x;
    if (abs_y > max_abs) max_abs = abs_y;
    if (abs_z > max_abs) max_abs = abs_z;
  }

  const double inv_n = 1.0 / static_cast<double>(n);

  c.x = (sum_x + comp_x) * inv_n;
  c.y = (sum_y + comp_y) * inv_n;
  c.z = (sum_z + comp_z) * inv_n;

  if (n > 1) {
    // Clamp tiny values to zero to avoid floating-point artifacts from cancellation
    constexpr double ABS_EPSILON = 1e-6;
    constexpr double REL_EPSILON = 1e-8;
    const double rel_eps         = REL_EPSILON * max_abs;
    const double eps             = (rel_eps > ABS_EPSILON) ? rel_eps : ABS_EPSILON;
    if (std::abs(c.x) < eps) c.x = 0.0;
    if (std::abs(c.y) < eps) c.y = 0.0;
    if (std::abs(c.z) < eps) c.z = 0.0;
  }
  return c;
}

} // namespace detail

struct AtomRec {
  AtomType type;
  // Technically, atoms are stores as T*, but I'm not sure upstream will work
  // with a nullable object. So we use T& as a mandatory never-null.
  std::reference_wrapper<const RDKit::Atom> atom;
};

struct RingRec {
  std::vector<std::reference_wrapper<const RDKit::Atom>> atoms;
  bool aromatic;

  RDGeom::Point3D center(const RDKit::Conformer &conf) const {
    return detail::centroid_from_atoms(atoms, conf);
  }

  RDGeom::Point3D normal(const RDKit::Conformer &conf) const {
    RDGeom::Point3D nrm{0.0, 0.0, 0.0};
    if (atoms.size() < 3) return nrm;

    const auto &a0 = atoms[0].get();
    const auto &a1 = atoms[1].get();
    const auto &a2 = atoms[2].get();

    auto p0 = conf.getAtomPos(a0.getIdx());
    auto p1 = conf.getAtomPos(a1.getIdx());
    auto p2 = conf.getAtomPos(a2.getIdx());
    auto v1 = p1 - p0;
    auto v2 = p2 - p0;
    nrm     = v1.crossProduct(v2);

    const double len_sq = nrm.lengthSq();
    constexpr double NORMAL_EPSILON = 1e-12;
    if (len_sq <= NORMAL_EPSILON) return RDGeom::Point3D{0.0, 0.0, 0.0};
    return nrm / std::sqrt(len_sq);
  }
};

struct GroupRec {
  AtomType a_type;
  FeatureGroup type;
  std::vector<std::reference_wrapper<const RDKit::Atom>> atoms;

  RDGeom::Point3D center(const RDKit::Conformer &conf) const {
    return detail::centroid_from_atoms(atoms, conf);
  }
};

template <typename T>
struct KindOf;

template <>
struct KindOf<AtomRec> {
  static constexpr Kind value = Kind::Atom;
};
template <>
struct KindOf<RingRec> {
  static constexpr Kind value = Kind::Ring;
};
template <>
struct KindOf<GroupRec> {
  static constexpr Kind value = Kind::Group;
};

template <Kind K>
struct RecordTypeFor;
template <>
struct RecordTypeFor<Kind::Atom> {
  using type = AtomRec;
};
template <>
struct RecordTypeFor<Kind::Ring> {
  using type = RingRec;
};
template <>
struct RecordTypeFor<Kind::Group> {
  using type = GroupRec;
};

} // namespace lahuta

#endif // LAHUTA_ENTITIES_RECORDS_HPP
