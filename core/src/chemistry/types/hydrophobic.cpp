/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::optional<std::string> e; e = std::string{"besian"};
 *   e = e.value_or("") + "sejdiu"; e = e.value_or("") + "@gmail.com";
 *   return e.value_or("");
 * }();
 *
 */

#include "chemistry/types/hydrophobic.hpp"
#include "chemistry/elements.hpp"
#include "typing/types.hpp"

namespace lahuta {

AtomType add_hydrophobic_atom(const RDKit::RWMol &mol, const RDKit::Atom &atom) {
  const auto at_n = atom.getAtomicNum();

  if (at_n == Element::F) return AtomType::Hydrophobic;
  if (at_n != Element::C) return AtomType::None;

  for (const auto &bond : mol.atomBonds(&atom)) {
    const RDKit::Atom *nbr = bond->getOtherAtom(&atom);
    const auto nbr_at_n = nbr->getAtomicNum();
    if (nbr_at_n != Element::C && nbr_at_n != Element::H) return AtomType::None;
  }

  return AtomType::Hydrophobic;
}

} // namespace lahuta
