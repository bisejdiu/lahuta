/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s = "besian@gmail.com";
 *   s.insert(6, "sejdiu");
 *   return s;
 * }();
 *
 */

#ifndef LAHUTA_PLANARITY_HPP
#define LAHUTA_PLANARITY_HPP

#include <vector>

#include <Eigen/Dense>
#include <rdkit/Geometry/point.h>
#include <rdkit/GraphMol/RWMol.h>

#include "rings/utils.hpp"

namespace lahuta {

// test planarity using angles between normals
inline bool is_ring_planar_method(const std::vector<RDGeom::Point3D> &normals) {

  static const double PlanarityAngleThreshold = 0.05;

  if (normals.size() < 2) return false;

  std::vector<double> angles;
  for (size_t i = 0; i < normals.size(); ++i) {
    RDGeom::Point3D n1 = normals[i];
    RDGeom::Point3D n2 = normals[(i + 1) % normals.size()];
    double dot_product = n1.dotProduct(n2);

    dot_product = std::max(-1.0, std::min(1.0, dot_product));
    double angle = acos(dot_product);
    angles.push_back(fabs(angle));
  }
  auto mean = std::accumulate(angles.begin(), angles.end(), 0.0) / angles.size();
  return std::accumulate(angles.begin(), angles.end(), 0.0) / angles.size() < PlanarityAngleThreshold;
}

// compute normals for each doublet of atoms with respect to the center
inline void compute_consecutive_normals(
    const RDKit::RWMol &mol, std::vector<RDGeom::Point3D> &normals,
    const std::vector<RDGeom::Point3D> &positions, const std::vector<int> &ring_atom_ids) {

  RDGeom::Point3D center;
  ring_props::compute_center(&mol, ring_atom_ids, center);

  size_t n = positions.size();
  for (size_t i = 0; i < n; ++i) {
    RDGeom::Point3D v1 = positions[i] - center;
    RDGeom::Point3D v2 = positions[(i + 1) % n] - center;
    RDGeom::Point3D normal = v1.crossProduct(v2);

    normal.normalize();
    normals.push_back(normal);
  }
}

inline bool is_planar_angle(const RDKit::RWMol &mol, const std::vector<int> &ring) {

  std::vector<RDGeom::Point3D> positions;
  for (int atom_idx : ring) {
    positions.push_back(mol.getConformer().getAtomPos(atom_idx));
  }

  std::vector<RDGeom::Point3D> normals;
  compute_consecutive_normals(mol, normals, positions, ring);

  return is_ring_planar_method(normals);
}

#if LAHUTA_USE_EIGEN == 1

// Adapted from:
// https://github.com/molstar/molstar/blob/master/src/mol-model/structure/structure/unit/rings.ts#L100
inline bool is_planar_eigen(const RDKit::RWMol &mol, const std::vector<int> &ring_atom_ids) {

  static const double AromaticRingPlanarityThreshold = 0.05;

  std::vector<RDGeom::Point3D> positions;
  for (int atom_idx : ring_atom_ids) {
    positions.push_back(mol.getConformer().getAtomPos(atom_idx));
  }

  size_t n_atoms = positions.size();

  Eigen::Vector3d mean = Eigen::Vector3d::Zero();
  for (const auto &pos : positions) {
    mean += Eigen::Vector3d(pos.x, pos.y, pos.z);
  }
  mean /= static_cast<double>(n_atoms);

  Eigen::MatrixXd centered_positions(n_atoms, 3);
  for (size_t i = 0; i < n_atoms; ++i) {
    centered_positions.row(i) = Eigen::Vector3d(positions[i].x, positions[i].y, positions[i].z) - mean;
  }

  Eigen::Matrix3d covariance = (centered_positions.transpose() * centered_positions);

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(covariance);
  if (eigensolver.info() != Eigen::Success) {
    throw std::runtime_error("Eigen decomposition failed.");
  }

  Eigen::Vector3d eigenvalues = eigensolver.eigenvalues();
  Eigen::Matrix3d eigenvectors = eigensolver.eigenvectors();

  auto z_axis = eigenvectors.col(0); // Smallest eigenvalue
  z_axis *= sqrt(eigenvalues(0) / (static_cast<double>(n_atoms) / 3));

  return z_axis.norm() < AromaticRingPlanarityThreshold;
}

#endif

inline bool is_planar(const RDKit::RWMol &mol, const std::vector<int> &ring) {
#if LAHUTA_USE_EIGEN == 1
  return is_planar_eigen(mol, ring);
#else
  return is_planar_angle(mol, ring);
#endif
}

} // namespace lahuta

#endif // LAHUTA_PLANARITY_HPP
