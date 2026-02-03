/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
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

void update_explicit_h_count(RDKit::Atom *atom_a, RDKit::Atom *atom_b) {
  int is_a_h = atom_a->getAtomicNum() == 1;
  int is_b_h = atom_b->getAtomicNum() == 1;

  // If exactly one atom is hydrogen, increment explicit H count on the non-hydrogen atom
  if (is_a_h ^ is_b_h) {
    auto non_h_atom = atom_a->getAtomicNum() == 1 ? atom_b : atom_a;
    non_h_atom->setNumExplicitHs(non_h_atom->getNumExplicitHs() + 1);
  }
}

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
      auto a = target.getAtomWithIdx(b_idx);
      auto b = target.getAtomWithIdx(e_idx);

      update_explicit_h_count(a, b); // Update explicit hydrogen count if forming H-X bond

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
