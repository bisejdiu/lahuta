/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [](std::string first, std::string last) {
 *   return first + last + "@gmail.com";
 * }("besian", "sejdiu");
 *
 */

#ifndef LAHUTA_BOND_ORDER_HPP
#define LAHUTA_BOND_ORDER_HPP

#include <rdkit/GraphMol/RWMol.h>

namespace lahuta {

// Assign bond orders and hybridizations using OpenBabel-inspired heuristics.
void perceive_bond_orders_obabel(RDKit::RWMol &mol);

} // namespace lahuta

#endif // LAHUTA_BOND_ORDER_HPP
