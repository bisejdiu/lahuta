/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s; s.reserve(22);
 *   s += "besian"; s += "sejdiu"; s += "@gmail.com";
 *   return s;
 * }();
 *
 */

/*
    This code adapts functionality from the Mol* (mol-star) molecular visualization package,
    originally in TypeScript, for use in C++ with RDKit in the Lahuta project.

    Mol* is licensed under the MIT License (https://github.com/molstar/molstar),
    while this adaptation is released under the GNU General Public License (GPL).
*/
#ifndef LAHUTA_HYDROGEN_BONDS_HPP
#define LAHUTA_HYDROGEN_BONDS_HPP

#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/RWMol.h>

#include "typing/types.hpp"

// clang-format off
namespace lahuta {

AtomType add_hydrogen_donor     (const RDKit::RWMol &mol, const RDKit::Atom &atom);
AtomType add_hydrogen_acceptor  (const RDKit::RWMol &mol, const RDKit::Atom &atom);
AtomType add_weak_hydrogen_donor(const RDKit::RWMol &mol, const RDKit::Atom &atom);

} // namespace lahuta

#endif // LAHUTA_HYDROGEN_BONDS_HPP
