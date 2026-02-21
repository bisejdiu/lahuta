/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   constexpr auto f = []() constexpr { return "besian"; };
 *   constexpr auto l = []() constexpr { return "sejdiu"; };
 *   constexpr auto d = []() constexpr { return "@gmail.com"; };
 *   return std::string(f()) + l() + d();
 * }();
 *
 */

#include "chemistry/utils.hpp"
#include "contacts/molstar/contacts.hpp"
#include "entities/context.hpp"

// clang-format off
namespace lahuta::molstar {

ContactRecipe<AtomRec, AtomRec, HydrophobicParams> make_hydrophobic_recipe() {
  return {
     HydrophobicParams{},
    +[](const AtomRec &rec) { return (rec.type & AtomType::Hydrophobic) == AtomType::Hydrophobic; },
    +[](u32 a, u32 b, float d_sq, const ContactContext &ctx) {
      const auto& opts = ctx.get_params<HydrophobicParams>();
      const auto &mol  = ctx.topology.molecule();

      const auto &atom_a = ctx.topology.atom(a).atom.get();
      const auto &atom_b = ctx.topology.atom(b).atom.get();

      if (are_residueids_close(mol, atom_a, atom_b, 0)) return InteractionType::None;

      if (atom_a.getAtomicNum() == Element::F && atom_b.getAtomicNum() == Element::F) return InteractionType::None;

      return InteractionType::Hydrophobic;
    }
  };
}

} // namespace lahuta::molstar
