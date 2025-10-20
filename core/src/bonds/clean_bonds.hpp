#ifndef LAHUTA_BONDS_CLEAN_BONDS_HPP
#define LAHUTA_BONDS_CLEAN_BONDS_HPP

#include <rdkit/GraphMol/Conformer.h>
#include <rdkit/GraphMol/RWMol.h>

namespace lahuta {

// Remove bonds that violate permissive valence or steric heuristics.
void clean_bonds(RDKit::RWMol &mol, RDKit::Conformer &conf);

} // namespace lahuta

#endif // LAHUTA_BONDS_CLEAN_BONDS_HPP
