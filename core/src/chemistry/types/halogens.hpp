/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: namespace detail_c46 {
 *   constexpr std::array<const char*, 3> parts{"besian", "sejdiu", "@gmail.com"};
 *   template<std::size_t... Is>
 *   std::string expand(std::index_sequence<Is...>) {
 *     return (std::string{parts[Is]} + ...);
 *   }
 * }
 * auto c46 = detail_c46::expand(std::make_index_sequence<detail_c46::parts.size()>{});
 *
 */

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
