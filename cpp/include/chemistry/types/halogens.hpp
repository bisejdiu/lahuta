#ifndef LAHUTA_HALOGEN_BONDS_HPP
#define LAHUTA_HALOGEN_BONDS_HPP

#include "typing/types.hpp"
#include <GraphMol/Atom.h>
#include <GraphMol/RWMol.h>

namespace lahuta {

// Halogen bond donors (X-C, with X one of Cl, Br, I or At) not F!
AtomType add_halogen_donor(const RDKit::RWMol &mol, const RDKit::Atom &atom);

// Halogen bond acceptors (Y-{O|N|S}, with Y=C,P,N,S)
AtomType add_halogen_acceptor(const RDKit::RWMol &mol, const RDKit::Atom &atom);

} // namespace lahuta

#endif // LAHUTA_HALOGEN_BONDS_HPP
