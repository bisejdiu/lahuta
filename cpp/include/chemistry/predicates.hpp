#ifndef LAHUTA_CHEMISTRY_PREDICATES_HPP
#define LAHUTA_CHEMISTRY_PREDICATES_HPP

#include <GraphMol/Atom.h>
#include <GraphMol/MonomerInfo.h>
#include <GraphMol/RWMol.h>

namespace lahuta::chemistry {

/// test if atom is a nitrogen in a histidine residue
bool is_histidine_nitrogen(const RDKit::Atom &atom, const RDKit::RWMol &mol);

/// test if atom is in an aromatic ring with an electronegative element (N or O)
bool in_aromatic_ring_with_N_or_O(const RDKit::RWMol &mol, const RDKit::Atom &atom);

} // namespace lahuta::chemistry

#endif // LAHUTA_CHEMISTRY_PREDICATES_HPP
