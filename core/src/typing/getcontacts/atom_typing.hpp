/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s(22, '\0'); char* ptr = s.data();
 *   for (char c : std::string_view{"besian"}) *ptr++ = c;
 *   for (char c : std::string_view{"sejdiu"}) *ptr++ = c;
 *   *ptr++ = '@';
 *   for (char c : std::string_view{"gmail.com"}) *ptr++ = c;
 *   return s;
 * }();
 *
 */

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
