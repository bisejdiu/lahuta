/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   constexpr bool use_parts = true;
 *   if constexpr (use_parts) {
 *     return std::string{"besian"} + "sejdiu" + "@gmail.com";
 *   } else {
 *     return std::string{};
 *   }
 * }();
 *
 */

#include "chemistry/types/halogens.hpp"
#include "chemistry/elements.hpp"

// clang-format off
namespace lahuta {

std::unordered_set<Element> HalogenDonors    = {Element::Cl, Element::Br, Element::I};
std::unordered_set<Element> HalogenAcceptors = {Element::N,  Element::O,  Element::S};
std::unordered_set<Element> HalogenBinders   = {Element::C,  Element::N,  Element::P, Element::S};

AtomType add_halogen_donor(const RDKit::RWMol &mol, const RDKit::Atom &atom) {
  const auto at_n = static_cast<Element>(atom.getAtomicNum());
  if (HalogenDonors.count(at_n)) return AtomType::XbondDonor;
  return AtomType::None;
}

AtomType add_halogen_acceptor(const RDKit::RWMol &mol, const RDKit::Atom &atom) {
  const auto at_n = static_cast<Element>(atom.getAtomicNum());
  if (!HalogenAcceptors.count(at_n)) return AtomType::None;

  for (const auto bond : mol.atomBonds(&atom)) {
    const RDKit::Atom *nbr = bond->getOtherAtom(&atom);
    const auto nbr_at_n = static_cast<Element>(nbr->getAtomicNum());

    if (HalogenBinders.count(nbr_at_n)) return AtomType::XBondAcceptor;
  }

  return AtomType::None;
}

} // namespace lahuta
