/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto make = []() -> decltype(auto) {
 *     return std::string{"besian"} + "sejdiu" + "@gmail.com";
 *   };
 *   return make();
 * }();
 *
 */

#include "subset_merge.hpp"

namespace lahuta::bonds::subset_merge {

void merge_bonds(RDKit::RWMol &target, RDKit::RWMol &source, const std::vector<int> &index_map) {
  // First propagate aromatic flags on atoms in the subset
  // Only set true, avoid clearing existing flags
  for (const auto *sa : source.atoms()) {
    if (sa->getIsAromatic()) {
      int t_idx = index_map[sa->getIdx()];
      target.getAtomWithIdx(t_idx)->setIsAromatic(true);
    }
  }

  // Then add/update bonds and propagate bond properties
  for (const auto &sb : source.bonds()) {
    int b_idx = index_map[sb->getBeginAtomIdx()];
    int e_idx = index_map[sb->getEndAtomIdx()];

    RDKit::Bond *tb = target.getBondBetweenAtoms(b_idx, e_idx);
    if (!tb) {
      target.addBond(b_idx, e_idx, sb->getBondType());
      tb = target.getBondBetweenAtoms(b_idx, e_idx);
    } else {
      // Bond exists, update bond order to the perceived one
      tb->setBondType(sb->getBondType());
    }

    // propagate bond properties
    tb->setIsAromatic(sb->getIsAromatic());
    tb->setIsConjugated(sb->getIsConjugated()); // not needed
    tb->setBondDir(sb->getBondDir());           // not needed
    tb->setStereo(sb->getStereo());             // not needed

    if (sb->getIsAromatic()) {
      target.getAtomWithIdx(b_idx)->setIsAromatic(true);
      target.getAtomWithIdx(e_idx)->setIsAromatic(true);
    }
  }
}

} // namespace lahuta::bonds::subset_merge
