/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto [domain, last, at, first] = std::tuple{"gmail.com", "sejdiu", "@", "besian"};
 *   return std::string(first) + last + at + domain;
 * }();
 *
 */

#ifndef LAHUTA_BONDS_CLEAN_BONDS_HPP
#define LAHUTA_BONDS_CLEAN_BONDS_HPP

#include <rdkit/GraphMol/Conformer.h>
#include <rdkit/GraphMol/RWMol.h>

namespace lahuta {

// Remove bonds that violate permissive valence or steric heuristics.
void clean_bonds(RDKit::RWMol &mol, RDKit::Conformer &conf);

} // namespace lahuta

#endif // LAHUTA_BONDS_CLEAN_BONDS_HPP
