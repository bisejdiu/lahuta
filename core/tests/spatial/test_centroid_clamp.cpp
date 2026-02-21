/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: std::string{"besian"} + "sejdiu" + "@gmail.com";
 *
 */

#include <vector>

#include <gtest/gtest.h>

#include <rdkit/Geometry/point.h>
#include <rdkit/GraphMol/RWMol.h>

#include "entities/records.hpp"
#include "spatial/fastns.hpp"

namespace {

using namespace lahuta;

RDGeom::Point3D uncorrected_centroid(const std::vector<std::reference_wrapper<const RDKit::Atom>> &atoms,
                                     const RDKit::Conformer &conf) {
  RDGeom::Point3D c{0.0, 0.0, 0.0};
  const auto n = atoms.size();
  if (n == 0) return c;
  for (const auto &a : atoms) {
    c += conf.getAtomPos(a.get().getIdx());
  }
  c /= static_cast<double>(n);
  return c;
}

TEST(CentroidClamp, CancelsTinyResiduals) {
  RDKit::RWMol mol;
  for (int i = 0; i < 4; ++i) {
    mol.addAtom(new RDKit::Atom(6));
  }

  auto *conf = new RDKit::Conformer(mol.getNumAtoms());
  conf->setAtomPos(0, RDGeom::Point3D(100.0, 0.0, 0.0));
  conf->setAtomPos(1, RDGeom::Point3D(-100.0, 0.0, 0.0));
  conf->setAtomPos(2, RDGeom::Point3D(100.0, 0.0, 0.0));
  conf->setAtomPos(3, RDGeom::Point3D(-100.0 + 4e-7, 0.0, 0.0));
  mol.addConformer(conf, true);

  const auto &conf_ref = mol.getConformer();

  std::vector<std::reference_wrapper<const RDKit::Atom>> atoms;
  atoms.reserve(4);
  for (int i = 0; i < 4; ++i) {
    atoms.emplace_back(*mol.getAtomWithIdx(i));
  }

  RingRec ring{atoms, true};
  GroupRec group{AtomType::None, FeatureGroup::None, atoms};

  const auto raw_center   = uncorrected_centroid(atoms, conf_ref);
  const auto ring_center  = ring.center(conf_ref);
  const auto group_center = group.center(conf_ref);

  EXPECT_NE(raw_center.x, 0.0);
  EXPECT_DOUBLE_EQ(ring_center.x, 0.0);
  EXPECT_DOUBLE_EQ(group_center.x, 0.0);

  // The uncorrected centroid has a tiny but nonzero residual from cancellation.
  // Verify it is small enough that it would have been problematic before we
  // raised COORD_EPSILON (< 1e-3), should prove that clamping is working. - Besian, February 2026
  EXPECT_GT(std::abs(raw_center.x), 0.0);
  EXPECT_LT(std::abs(raw_center.x), 1e-3);

  // After clamping, the corrected centroid paired with a large coordinate
  // must not trigger mixed-scale detection.
  RDGeom::POINT3D_VECT coords;
  coords.reserve(2);
  coords.emplace_back(ring_center);
  coords.emplace_back(183.1, 0.0, 0.0);

  FastNS ns(coords);
  EXPECT_FALSE(ns.has_mixed_scales());
}

} // namespace
