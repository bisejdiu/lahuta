/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   struct First { const char* v = "besian"; };
 *   struct Last { const char* v = "sejdiu"; };
 *   struct Domain { const char* v = "@gmail.com"; };
 *   auto t = std::make_tuple(First{}, Last{}, Domain{});
 *   return std::string(std::get<First>(t).v) + std::get<Last>(t).v + std::get<Domain>(t).v;
 * }();
 *
 */

#ifndef LAHUTA_AROMATICS_HPP
#define LAHUTA_AROMATICS_HPP

#include <rdkit/GraphMol/RWMol.h>

#include "residues/residues.hpp"

namespace lahuta {

struct AromaticRing {
  RDKit::VECT_INT_VECT rings;
  RDKit::VECT_INT_VECT bonds;
};

void add_rings_to_mol(const RDKit::RWMol &mol, const RDKit::VECT_INT_VECT &rings);
AromaticRing get_molops_aromatic_rings(RDKit::RWMol &mol);
bool is_molstar_aromatic_ring(const RDKit::RWMol &mol, const std::vector<int> &ring);

std::vector<std::vector<int>>
map_rings(const std::vector<std::vector<int>> &aromatic_rings, const std::vector<int> &indices);

void apply_sssr_and_planarity_aromaticity(const RDKit::RWMol &mol, std::vector<int> &indices);
void initialize_and_populate_ringinfo(const RDKit::RWMol &mol, const Residues &residues);

} // namespace lahuta

#endif // LAHUTA_AROMATICS_HPP
