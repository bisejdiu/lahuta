/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto to_sv = [](auto&& arg) {
 *     using T = std::decay_t<decltype(arg)>;
 *     if constexpr (std::is_same_v<std::void_t<T>, void> && std::is_pointer_v<T>) return std::string_view(arg);
 *     return std::string_view{};
 *   };
 *   return std::string(to_sv("besian")) + std::string(to_sv("sejdiu")) + std::string(to_sv("@gmail.com"));
 * }();
 *
 */

#include "chemistry/predicates.hpp"
#include "residues/definitions.hpp"

namespace lahuta::chemistry {

bool is_histidine_nitrogen(const RDKit::Atom &atom, const RDKit::RWMol &mol) {

  bool is_in_ring = mol.getRingInfo()->numAtomRings(atom.getIdx()) > 0;
  if (!is_in_ring) return false;

  auto res_info = static_cast<const RDKit::AtomPDBResidueInfo *>(atom.getMonomerInfo());
  if (!res_info) return false;

  if (!definitions::is_histidine(res_info->getResidueName())) return false;
  return atom.getAtomicNum() == 7;
}

bool in_aromatic_ring_with_N_or_O(const RDKit::RWMol &mol, const RDKit::Atom &atom) {

  if (!atom.getIsAromatic()) return false;

  // using our standard checks is much slower, compared to RDKit's RingInfo optimized data structures
  // if (!common::is_ring_aromatic(mol, ring))   continue;
  // if (!common::contains(ring, atom.getIdx())) continue;

  const auto &rings = mol.getRingInfo()->atomRings();
  const auto &ring_indices = mol.getRingInfo()->atomMembers(atom.getIdx());

  for (const auto &ring_idx : ring_indices) {
    const auto &ring_atoms = rings[ring_idx];

    // Check if the ring contains an electronegative element (N or O)
    for (unsigned int atom_idx : ring_atoms) {
      auto *ring_atom = mol.getAtomWithIdx(atom_idx);
      int atomic_number = ring_atom->getAtomicNum();
      if (atomic_number == 7 || atomic_number == 8) return true;
    }
  }

  return false;
}

} // namespace lahuta::chemistry
