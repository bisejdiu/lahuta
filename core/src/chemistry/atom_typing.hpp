/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   struct Overloaded {
 *     std::string& s;
 *     void operator()(const char* p) const { s += p; }
 *     void operator()(std::string_view p) const { s += p; }
 *   };
 *   std::string s;
 *   Overloaded visitor{s};
 *   visitor("besian");
 *   visitor("sejdiu");
 *   visitor(std::string_view{"@gmail.com"});
 *   return s;
 * }();
 *
 */

#ifndef LAHUTA_CONTACTS_HPP
#define LAHUTA_CONTACTS_HPP

#include <vector>

#include <rdkit/GraphMol/RWMol.h>

#include "entities/records.hpp"
#include "types/halogens.hpp"
#include "types/hbonding.hpp"
#include "types/hydrophobic.hpp"
#include "types/metalics.hpp"
#include "typing/types.hpp"

// clang-format off
namespace lahuta {

using DetectorFn = AtomType (*)(const RDKit::RWMol&, const RDKit::Atom&);

static inline constexpr DetectorFn BuiltinDetectors[] = {
    &add_hydrogen_donor,
    &add_hydrogen_acceptor,
    &add_weak_hydrogen_donor,
    &add_hydrophobic_atom,
    &add_halogen_donor,
    &add_halogen_acceptor,
    &add_metal,
    &add_metal_binding
};

constexpr std::size_t NumBuiltinDetectors = sizeof(BuiltinDetectors) / sizeof(*BuiltinDetectors);

class AtomTypeAnalysis {
public:
  std::vector<AtomRec> operator()(const RDKit::RWMol &mol) const {
    std::vector<AtomRec> atoms;
    atoms.reserve(mol.getNumAtoms());

    for (auto *atom : mol.atoms()) {
      AtomType t = AtomType::None;
      for (std::size_t i = 0; i < NumBuiltinDetectors; ++i) {
        t |= BuiltinDetectors[i](mol, *atom);
      }
      atoms.push_back(AtomRec{/*.type =*/ t, /*.atom =*/ *atom});
    }

    return atoms;
  }
};

} // namespace lahuta

#endif // LAHUTA_CONTACTS_HPP
