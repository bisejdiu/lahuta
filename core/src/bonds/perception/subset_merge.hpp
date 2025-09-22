#pragma once

#include <vector>

#include <rdkit/GraphMol/RWMol.h>

namespace lahuta::bonds::subset_merge {

void merge_bonds(RDKit::RWMol &target, RDKit::RWMol &source, const std::vector<int> &index_map);
void update_explicit_h_count(RDKit::Atom *atom_a, RDKit::Atom *atom_b);

} // namespace lahuta::bonds::subset_merge
