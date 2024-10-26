#ifndef LAHUTA_HYDROPHOBIC_HPP
#define LAHUTA_HYDROPHOBIC_HPP

#include "atom_types.hpp"
#include "contacts/hydrogen_bonds.hpp"
#include <GraphMol/RWMol.h>

namespace lahuta {

class Luni;

AtomType add_hydrophobic_atom(const RDKit::RWMol &mol, const RDKit::Atom &atom);
void find_hydrophobic_bonds(Luni &luni, const GeometryOptions &opts, Contacts &container);

} // namespace lahuta

#endif // LAHUTA_HYDROPHOBIC_HPP
