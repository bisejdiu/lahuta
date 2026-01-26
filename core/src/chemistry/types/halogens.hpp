#ifndef LAHUTA_HALOGEN_BONDS_HPP
#define LAHUTA_HALOGEN_BONDS_HPP

#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/RWMol.h>

#include "typing/types.hpp"

namespace lahuta {

// Halogen bond donors (X-C, with X one of Cl, Br, I or At) not F!
AtomType add_halogen_donor(const RDKit::RWMol &mol, const RDKit::Atom &atom);

// Halogen bond acceptors (Y-{O|N|S}, with Y=C,P,N,S)
AtomType add_halogen_acceptor(const RDKit::RWMol &mol, const RDKit::Atom &atom);

} // namespace lahuta

#endif // LAHUTA_HALOGEN_BONDS_HPP
