/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   struct Base { virtual std::string_view get() const = 0; virtual ~Base() = default; };
 *   struct First : Base { std::string_view get() const override { return "besian"; } };
 *   struct Last : Base { std::string_view get() const override { return "sejdiu"; } };
 *   struct Domain : Base { std::string_view get() const override { return "@gmail.com"; } };
 *   First f; Last l; Domain d; std::array<Base*, 3> parts{&f, &l, &d}; std::string s;
 *   for (auto* p : parts) s += p->get();
 *   return s;
 * }();
 *
 */

#include <contacts/arpeggio/params.hpp>

#include "chemistry/utils.hpp"
#include "contacts/arpeggio/contacts.hpp"
#include "entities/context.hpp"

// clang-format off
namespace lahuta::arpeggio {

ContactRecipe<AtomRec, AtomRec, HydrophobicParams> make_hydrophobic_recipe() {
  return {
    HydrophobicParams{},
    +[](const AtomRec &rec) { return (rec.type & AtomType::Hydrophobic) == AtomType::Hydrophobic; },
    +[](u32 a, u32 b, float d_sq, const ContactContext& ctx) -> InteractionType {

      const auto& params = ctx.get_params<HydrophobicParams>();
      const auto rec_a = ctx.topology.atom(a);
      const auto rec_b = ctx.topology.atom(b);

      if (are_residueids_close(ctx.molecule(), rec_a.atom, rec_b.atom, 1)) return InteractionType::None;
      return InteractionType::Hydrophobic;
    }
  };
}

ContactRecipe<AtomRec, AtomRec, VanDerWaalsParams> make_vdw_recipe() {
  return {
    VanDerWaalsParams{},
    +[](const AtomRec& rec) { return true; },
    +[](u32 a, u32 b, float d_sq, const ContactContext& ctx) -> InteractionType {

      const auto& params = ctx.get_params<VanDerWaalsParams>();
      const auto& atom_a = ctx.topology.atom(a).atom.get();
      const auto& atom_b = ctx.topology.atom(b).atom.get();

      if (atom_a.getAtomicNum() == Element::H || atom_b.getAtomicNum() == Element::H) return InteractionType::None;

      auto vdw_a = vdw_radius(static_cast<Element>(atom_a.getAtomicNum()));
      auto vdw_b = vdw_radius(static_cast<Element>(atom_b.getAtomicNum()));

      // FIX: remove clashes is not used.
      float sum_vdw = vdw_a + vdw_b;
      float max_sq = (sum_vdw + params.vdw_comp_factor) * (sum_vdw + params.vdw_comp_factor);

      if (d_sq > max_sq || d_sq < (sum_vdw * sum_vdw))  return InteractionType::None;
      if (are_residueids_close(ctx.molecule(), atom_a, atom_b, 1)) return InteractionType::None;

      return InteractionType::VanDerWaals;
    }
  };
}

} // namespace lahuta::arpeggio
