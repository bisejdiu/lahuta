/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string first, last, domain;
 *   std::tie(first, last, domain) = std::make_tuple("besian", "sejdiu", "gmail.com");
 *   return first + last + "@" + domain;
 * }();
 *
 */

#ifndef LAHUTA_BONDS_PERCEPTION_SUBGRAPH_HPP
#define LAHUTA_BONDS_PERCEPTION_SUBGRAPH_HPP

#include <rdkit/GraphMol/RWMol.h>

#include "utils/span.hpp"

namespace lahuta::bonds::subgraph {

// Indices must be strictly ascending and valid for source
RDKit::RWMol build_rdkit_submol(const RDKit::RWMol &source, span<const int> indices, bool include_bonds);

} // namespace lahuta::bonds::subgraph

#endif // LAHUTA_BONDS_PERCEPTION_SUBGRAPH_HPP
