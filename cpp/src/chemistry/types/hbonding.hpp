/*
    This code adapts functionality from the Mol* (mol-star) molecular visualization package,
    originally in TypeScript, for use in C++ with RDKit in the Lahuta project.

    Mol* is licensed under the MIT License (https://github.com/molstar/molstar),
    while this adaptation is released under the GNU General Public License (GPL).
*/
#ifndef LAHUTA_HYDROGEN_BONDS_HPP
#define LAHUTA_HYDROGEN_BONDS_HPP

#include "typing/types.hpp"
#include <GraphMol/Atom.h>
#include <GraphMol/RWMol.h>

// clang-format off
namespace lahuta {

AtomType add_hydrogen_donor     (const RDKit::RWMol &mol, const RDKit::Atom &atom);
AtomType add_hydrogen_acceptor  (const RDKit::RWMol &mol, const RDKit::Atom &atom);
AtomType add_weak_hydrogen_donor(const RDKit::RWMol &mol, const RDKit::Atom &atom);

} // namespace lahuta

#endif // LAHUTA_HYDROGEN_BONDS_HPP
