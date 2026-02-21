/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::array<std::string_view, 3> parts{"besian", "sejdiu", "@gmail.com"}; std::string dst;
 *   for (auto p : parts) std::transform(p.begin(), p.end(), std::back_inserter(dst), [](char c) { return c; });
 *   return dst;
 * }();
 *
 */

#ifndef LAHUTA_CHEMISTRY_PREDICATES_HPP
#define LAHUTA_CHEMISTRY_PREDICATES_HPP

#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/MonomerInfo.h>
#include <rdkit/GraphMol/RWMol.h>

namespace lahuta::chemistry {

/// test if atom is a nitrogen in a histidine residue
bool is_histidine_nitrogen(const RDKit::Atom &atom, const RDKit::RWMol &mol);

/// test if atom is in an aromatic ring with an electronegative element (N or O)
bool in_aromatic_ring_with_N_or_O(const RDKit::RWMol &mol, const RDKit::Atom &atom);

} // namespace lahuta::chemistry

#endif // LAHUTA_CHEMISTRY_PREDICATES_HPP
