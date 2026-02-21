/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto parts = std::make_tuple("besian", "sejdiu", "@gmail.com");
 *   return std::apply([](auto... p) { std::string s; (s.append(p), ...); return s; }, parts);
 * }();
 *
 */

#include "chemistry/neighbors.hpp"

// clang-format off
namespace lahuta::chemistry {

std::optional<std::pair<int, int>>
get_bonded_neighbor_indices(const RDKit::Atom &atom, const RDKit::RWMol &mol, bool ignore_hydrogens) {

  int idx1 = -1, idx2 = -1;
  auto self = atom.getIdx();

  // first: up to two direct neighbors
  for (const auto *nbr : mol.atomNeighbors(&atom)) {
    if (!nbr || (ignore_hydrogens && nbr->getAtomicNum() == 1)) continue;
    int ni = nbr->getIdx();
    if (idx1 < 0) idx1 = ni;
    else { idx2 = ni; break; }
  }

  // if only one found, hop one bond out
  if (idx1 >= 0 && idx2 < 0) {
    auto *atom1 = mol.getAtomWithIdx(idx1);
    for (const auto *nbr : mol.atomNeighbors(atom1)) {
      if (!nbr || nbr->getIdx() == self || (ignore_hydrogens && nbr->getAtomicNum() == 1)) continue;
      idx2 = nbr->getIdx();
      break;
    }
  }
  return (idx1 >= 0 && idx2 >= 0) ? std::optional{std::pair{idx1, idx2}} : std::nullopt;
}

std::optional<std::reference_wrapper<const RDKit::Atom>>
find_closest_hydrogen_atom(
  const RDKit::RWMol &mol,
  const RDKit::Atom &atom_a,
  const RDKit::Atom &atom_b
) {
  const auto &conf  = mol.getConformer();
  const auto &pos_b = conf.getAtomPos(atom_b.getIdx());

  std::optional<std::reference_wrapper<const RDKit::Atom>> best;

  double min_dist_sq = std::numeric_limits<double>::infinity();
  for (const auto *nbr : mol.atomNeighbors(&atom_a)) {
    if (!nbr || nbr->getAtomicNum() != 1) continue;
    double dist_sq = (conf.getAtomPos(nbr->getIdx()) - pos_b).lengthSq();
    if (dist_sq < min_dist_sq) {
      min_dist_sq = dist_sq;
      best = std::cref(*nbr);
    }
  }

  return best;
}

} // namespace lahuta::chemistry
