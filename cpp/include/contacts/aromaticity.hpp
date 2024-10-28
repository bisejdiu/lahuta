#ifndef LAHUTA_AROMATICITY_HPP
#define LAHUTA_AROMATICITY_HPP

#include "GraphMol/RWMol.h"
#include "RDGeneral/types.h"
#include "contacts/features.hpp"
#include "ob/kekulize.h"
#include "residues.hpp"

namespace lahuta {

void add_rings_to_mol(const RDKit::RWMol &mol, const RDKit::VECT_INT_VECT &rings);

struct AromaticRing {
  RDKit::VECT_INT_VECT rings;
  RDKit::VECT_INT_VECT bonds;
};
AromaticRing get_molops_aromatic_rings(RDKit::RWMol &mol);

std::vector<const RDKit::Atom *> get_atoms(const RDKit::RWMol &mol, const std::vector<int> &indices);

// computes aromatic rings from filtered atoms (based on indices) to mol
FeatureVec add_aromatic_rings(const RDKit::RWMol &mol, const std::vector<int> &indices);

/// add aromatic rings
FeatureVec add_aromatic_rings(const RDKit::RWMol &mol, Residues &residues);

} // namespace lahuta

#endif // LAHUTA_AROMATICITY_HPP
