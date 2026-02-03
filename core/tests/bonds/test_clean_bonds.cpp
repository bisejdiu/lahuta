/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [](std::string first, std::string last) {
 *   return first + last + "@gmail.com";
 * }("besian", "sejdiu");
 *
 */

#include <cmath>
#include <memory>
#include <vector>

#include <gtest/gtest.h>

#include <rdkit/Geometry/point.h>
#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/Conformer.h>
#include <rdkit/GraphMol/RWMol.h>

#include "bonds/clean_bonds.hpp"

namespace lahuta {
namespace {

std::unique_ptr<RDKit::Conformer> make_conformer(const std::vector<RDGeom::Point3D> &coords) {
  auto conf = std::make_unique<RDKit::Conformer>(static_cast<unsigned>(coords.size()));
  conf->set3D(true);
  for (unsigned i = 0; i < coords.size(); ++i) {
    conf->setAtomPos(i, coords[i]);
  }
  return conf;
}

RDKit::RWMol build_molecule(const std::vector<unsigned> &atomic_numbers) {
  RDKit::RWMol mol;
  for (unsigned atomic_num : atomic_numbers) {
    mol.addAtom(new RDKit::Atom(atomic_num), false, true);
  }
  return mol;
}

} // namespace

TEST(CleanBondsTest, RemovesLongestBondWhenAboveValenceLimit) {
  // Carbon with five sigma bonds should drop the farthest one.
  auto mol = build_molecule({6, 1, 1, 1, 1, 1});
  for (unsigned idx = 1; idx <= 5; ++idx) {
    mol.addBond(0, idx, RDKit::Bond::BondType::SINGLE);
  }

  std::vector<RDGeom::Point3D> coords = {
      {0.0, 0.0, 0.0},
      {1.0, 0.0, 0.0},
      {-1.0, 0.0, 0.0},
      {0.0, 1.0, 0.0},
      {0.0, -1.0, 0.0},
      {0.0, 0.0, 3.0}, // furthest hydrogen, should be pruned
  };
  auto conf = make_conformer(coords);
  mol.addConformer(conf.release(), true);

  auto &active_conf = mol.getConformer();
  clean_bonds(mol, active_conf);

  EXPECT_EQ(mol.getAtomWithIdx(0)->getDegree(), 4u);
  EXPECT_NE(nullptr, mol.getBondBetweenAtoms(0, 1));
  EXPECT_EQ(nullptr, mol.getBondBetweenAtoms(0, 5));
}

TEST(CleanBondsTest, PrefersRemovingHydrogenHydrogenBonds) {
  // Over-coordinated hydrogen should discard the H-H contact first.
  auto mol = build_molecule({6, 1, 1});
  mol.addBond(0, 1, RDKit::Bond::BondType::SINGLE);
  mol.addBond(1, 2, RDKit::Bond::BondType::SINGLE);

  std::vector<RDGeom::Point3D> coords = {
      {0.0, 0.0, 0.0},
      {0.9, 0.0, 0.0},
      {1.8, 0.0, 0.0},
  };
  auto conf = make_conformer(coords);
  mol.addConformer(conf.release(), true);

  auto &active_conf = mol.getConformer();
  clean_bonds(mol, active_conf);

  EXPECT_NE(nullptr, mol.getBondBetweenAtoms(0, 1));
  EXPECT_EQ(nullptr, mol.getBondBetweenAtoms(1, 2));
  EXPECT_EQ(mol.getAtomWithIdx(1)->getDegree(), 1u);
}

TEST(CleanBondsTest, RemovesBondFormingAcuteAngle) {
  // Center atom keeps at most one of the nearly collinear bonds.
  auto mol = build_molecule({6, 7, 7, 8});
  mol.addBond(0, 1, RDKit::Bond::BondType::SINGLE);
  mol.addBond(0, 2, RDKit::Bond::BondType::SINGLE);
  mol.addBond(0, 3, RDKit::Bond::BondType::SINGLE);

  constexpr double deg_to_rad = 3.14159265358979323846 / 180.0;
  const double radians = 20.0 * deg_to_rad;
  std::vector<RDGeom::Point3D> coords = {
      {0.0, 0.0, 0.0},
      {1.0, 0.0, 0.0},
      {std::cos(radians), std::sin(radians), 0.0},
      {0.0, 0.0, 1.5},
  };
  auto conf = make_conformer(coords);
  mol.addConformer(conf.release(), true);

  auto &active_conf = mol.getConformer();
  clean_bonds(mol, active_conf);

  const unsigned remaining_degree = mol.getAtomWithIdx(0)->getDegree();
  EXPECT_EQ(remaining_degree, 2u);
  const auto *bond12 = mol.getBondBetweenAtoms(0, 1);
  const auto *bond13 = mol.getBondBetweenAtoms(0, 2);
  EXPECT_TRUE((bond12 == nullptr) ^ (bond13 == nullptr));
  EXPECT_NE(nullptr, mol.getBondBetweenAtoms(0, 3));
}

} // namespace lahuta
