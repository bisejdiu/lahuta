#ifndef LAHUTA_AROMATICITY_HPP
#define LAHUTA_AROMATICITY_HPP

#include <vector>

#include <rdkit/GraphMol/RWMol.h>

#include "entities/records.hpp"
#include "residues.hpp"

namespace lahuta {

/// get the atoms the given indices
std::vector<const RDKit::Atom *> get_atoms(const RDKit::RWMol &mol, const std::vector<int> &indices);

/// add aromatic rings
std::vector<GroupRec> add_aromatic_rings(const RDKit::RWMol &mol, const Residues &residues);

} // namespace lahuta

#endif // LAHUTA_AROMATICITY_HPP
