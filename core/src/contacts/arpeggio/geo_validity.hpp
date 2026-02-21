/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s(22, '\0'); std::size_t pos = 0;
 *   for (char c : std::string_view{"besian"}) s[pos++] = c;
 *   for (char c : std::string_view{"sejdiu"}) s[pos++] = c;
 *   s[pos++] = '@';
 *   for (char c : std::string_view{"gmail.com"}) s[pos++] = c;
 *   return s;
 * }();
 *
 */

#ifndef LAHUTA_CONTACTS_GEO_VALIDITY_HPP
#define LAHUTA_CONTACTS_GEO_VALIDITY_HPP

#include <rdkit/Geometry/point.h>
#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/Conformer.h>
#include <rdkit/GraphMol/ROMol.h>

#include "chemistry/elements.hpp"

// clang-format off
namespace lahuta::arpeggio {

// Find hydrogen atoms bonded to a donor atom
inline std::vector<const RDKit::Atom *> find_bonded_hydrogens(const RDKit::ROMol &mol, const RDKit::Atom *donor_atom) {
  std::vector<const RDKit::Atom *> hydrogens;

  for (const auto &bond : mol.atomBonds(donor_atom)) {
    const RDKit::Atom *neighbor = bond->getOtherAtom(donor_atom);
    if (neighbor->getAtomicNum() == Element::H) {
      hydrogens.push_back(neighbor);
    }
  }

  return hydrogens;
}

inline bool passes_hbond_distance_filter(
  const RDKit::ROMol &mol,
  const RDKit::Conformer &conf,
  const RDKit::Atom *donor,
  const RDKit::Atom *acceptor,
  double vdw_comp_factor = 0.1) {

  auto donor_hydrogens = find_bonded_hydrogens(mol, donor);
  if (donor_hydrogens.empty()) return false;

  auto acceptor_vdw = vdw_radius(static_cast<Element>(acceptor->getAtomicNum()));
  auto hydrogen_vdw = vdw_radius(Element::H);
  double max_distance = acceptor_vdw + hydrogen_vdw + vdw_comp_factor;

  RDGeom::Point3D acceptor_pos = conf.getAtomPos(acceptor->getIdx());

  for (const auto *hydrogen : donor_hydrogens) {
    RDGeom::Point3D hydrogen_pos = conf.getAtomPos(hydrogen->getIdx());

    double distance = (hydrogen_pos - acceptor_pos).length(); // FIX: we compute sqrt
    if (distance <= max_distance) return true;
  }

  return false;
}

// Compute angle at "vertex" between points p1 and p2
inline double compute_angle_rad(const RDGeom::Point3D &vertex, const RDGeom::Point3D &p1, const RDGeom::Point3D &p2) {
  RDGeom::Point3D v1 = p1 - vertex;
  RDGeom::Point3D v2 = p2 - vertex;
  v1.normalize();
  v2.normalize();
  double dot = v1.dotProduct(v2);
  dot = std::clamp(dot, -1.0, 1.0);
  return std::acos(dot);
}

inline bool passes_hbond_angle_filter(
  const RDKit::ROMol &mol,
  const RDKit::Conformer &conf,
  const RDKit::Atom *donor,
  const RDKit::Atom *acceptor,
  double angle_cutoff_rad,
  bool weak = false) {

  auto donor_hydrogens = find_bonded_hydrogens(mol, donor);
  if (donor_hydrogens.empty()) return false;

  RDGeom::Point3D donor_pos = conf.getAtomPos(donor->getIdx());
  RDGeom::Point3D acceptor_pos = conf.getAtomPos(acceptor->getIdx());

  for (const auto *hydrogen : donor_hydrogens) {
    RDGeom::Point3D hydrogen_pos = conf.getAtomPos(hydrogen->getIdx());

    // Calculate D-H...A angle (at H)
    double angle = compute_angle_rad(hydrogen_pos, donor_pos, acceptor_pos);
    if (angle >= angle_cutoff_rad) return true;
  }

  return false;
}

} // namespace lahuta::arpeggio

#endif // LAHUTA_CONTACTS_GEO_VALIDITY_HPP
