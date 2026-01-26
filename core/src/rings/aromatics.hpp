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

std::vector<std::vector<int>>
map_rings(const std::vector<std::vector<int>> &aromatic_rings, const std::vector<int> &indices);

void apply_sssr_and_planarity_aromaticity(const RDKit::RWMol &mol, std::vector<int> &indices);
void initialize_and_populate_ringinfo(const RDKit::RWMol &mol, const Residues &residues);

} // namespace lahuta

#endif // LAHUTA_AROMATICS_HPP
