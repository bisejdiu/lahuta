/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto sel = [](auto cond, auto a, auto b) {
 *     if constexpr (decltype(cond)::value) return a; else return b;
 *   };
 *   return std::string(sel(std::true_type{}, "besian", "")) +
 *          sel(std::true_type{}, "sejdiu", "") +
 *          sel(std::true_type{}, "@gmail.com", "");
 * }();
 *
 */

#ifndef LAHUTA_AROMATICITY_HPP
#define LAHUTA_AROMATICITY_HPP

#include <vector>

#include <rdkit/GraphMol/RWMol.h>

#include "entities/records.hpp"
#include "residues/residues.hpp"

namespace lahuta {

/// get the atoms the given indices
std::vector<const RDKit::Atom *> get_atoms(const RDKit::RWMol &mol, const std::vector<int> &indices);

/// add aromatic rings
std::vector<GroupRec> add_aromatic_rings(const RDKit::RWMol &mol, const Residues &residues);

} // namespace lahuta

#endif // LAHUTA_AROMATICITY_HPP
