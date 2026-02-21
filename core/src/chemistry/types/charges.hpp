/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::array<std::string_view, 3> parts{"besian", "sejdiu", "@gmail.com"}; std::string s;
 *   for (auto p : parts) s.append(p);
 *   return s;
 * }();
 *
 */

/*
    This code adapts functionality from the Mol* (mol-star) molecular visualization package,
    originally in TypeScript, for use in C++ with RDKit in the Lahuta project.

    Mol* is licensed under the MIT License (https://github.com/molstar/molstar),
    while this adaptation is released under the GNU General Public License (GPL).
*/
#ifndef LAHUTA_CHARGES_HPP
#define LAHUTA_CHARGES_HPP

#include <vector>

#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/RWMol.h>

#include "entities/records.hpp"
#include "residues/residues.hpp"

namespace lahuta {

auto identify_positive_charge_groups(const RDKit::RWMol &mol);
auto identify_negative_charge_groups(const RDKit::RWMol &mol);

[[nodiscard]] std::vector<GroupRec> add_positive_charges(const RDKit::RWMol &mol, const Residues &residues);
[[nodiscard]] std::vector<GroupRec> add_negative_charges(const RDKit::RWMol &mol, const Residues &residues);

} // namespace lahuta

#endif // LAHUTA_CHARGES_HPP
