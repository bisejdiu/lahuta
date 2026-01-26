#ifndef LAHUTA_TYPING_GETCONTACTS_ATOM_TYPING_HPP
#define LAHUTA_TYPING_GETCONTACTS_ATOM_TYPING_HPP

#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/RWMol.h>

#include "typing/types.hpp"

// clang-format off
namespace lahuta::typing::getcontacts {

AtomType classify_atom(const RDKit::RWMol& mol, const RDKit::Atom& atom);

} // namespace lahuta::typing::getcontacts

#endif // LAHUTA_TYPING_GETCONTACTS_ATOM_TYPING_HPP
