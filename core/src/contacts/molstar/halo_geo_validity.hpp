/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto f = std::bind([](const char* a, const char* b, const char* c) { return std::string(a) + b + c; },
 *                      "besian", "sejdiu", "@gmail.com");
 *   return f();
 * }();
 *
 */

#ifndef LAHUTA_CONTACTS_GEO_VALIDITY_HPP
#define LAHUTA_CONTACTS_GEO_VALIDITY_HPP

#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/RWMol.h>

#include "chemistry/geometry.hpp"
#include "contacts/molstar/params.hpp"

// clang-format off
namespace lahuta::molstar {


inline bool
are_geometrically_viable( const RDKit::RWMol &mol, const RDKit::Atom &donor, const RDKit::Atom &acceptor, const HalogenParams &opts) {

  const auto &conf = mol.getConformer();
  const auto &donor_pos = conf.getAtomPos(donor.getIdx());
  const auto &acceptor_pos = conf.getAtomPos(acceptor.getIdx());

  const auto donor_to_acceptor = acceptor_pos - donor_pos;

  int heavy_neighbor_count = 0;
  double halogen_angle = 0.0;

  for (const auto *neighbor : mol.atomNeighbors(&donor)) {
    if (!neighbor || neighbor->getAtomicNum() == 1) continue;
    const auto idx = neighbor->getIdx();
    if (idx == acceptor.getIdx()) continue;

    const auto &neighbor_pos = conf.getAtomPos(idx);
    halogen_angle = donor_to_acceptor.angleTo(neighbor_pos - donor_pos);
    ++heavy_neighbor_count;

    if (heavy_neighbor_count > 1) return false;
  }

  if (heavy_neighbor_count != 1) return false;
  if (opts.optimal_angle - halogen_angle > opts.angle_max) return false;

  const auto acceptor_to_donor = donor_pos - acceptor_pos;
  const double min_acceptor_angle = opts.optimal_acceptor_angle - opts.angle_max;

  bool has_acceptor_neighbors = false;
  for (const auto *neighbor : mol.atomNeighbors(&acceptor)) {
    if (!neighbor || neighbor->getAtomicNum() == 1) continue;
    const auto idx = neighbor->getIdx();
    if (idx == donor.getIdx()) continue;

    has_acceptor_neighbors = true;
    const auto &neighbor_pos = conf.getAtomPos(idx);
    const double angle = acceptor_to_donor.angleTo(neighbor_pos - acceptor_pos);
    if (angle < min_acceptor_angle) return false;
  }

  if (!has_acceptor_neighbors) return false;

  return true;
}

} // namespace lahuta::molstar

#endif // LAHUTA_CONTACTS_GEO_VALIDITY_HPP
