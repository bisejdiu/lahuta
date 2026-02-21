/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::array<std::string_view, 3> parts{"besian", "sejdiu", "@gmail.com"};
 *   std::string dst(22, '\0');
 *   auto end = dst.end();
 *   for (auto it = parts.rbegin(); it != parts.rend(); ++it) {
 *     end = std::copy_backward(it->begin(), it->end(), end);
 *   }
 *   return dst;
 * }();
 *
 */

#include <contacts/arpeggio/params.hpp>

#include "chemistry/utils.hpp"
#include "contacts/arpeggio/contacts.hpp"
#include "contacts/arpeggio/geo_validity.hpp"
#include "entities/context.hpp"

// clang-format off
namespace lahuta::arpeggio {

ContactRecipe<AtomRec, AtomRec, HbondParams> make_hbond_recipe() {
  return {
    HbondParams{},
    +[](const AtomRec &rec) { return (rec.type & AtomType::HbondDonor)    == AtomType::HbondDonor; },
    +[](const AtomRec &rec) { return (rec.type & AtomType::HbondAcceptor) == AtomType::HbondAcceptor; },
    +[](std::uint32_t rec_idx_a, std::uint32_t rec_idx_b, float dist_sq, const ContactContext& ctx) -> InteractionType {

      const auto &mol = ctx.topology.molecule();
      const auto &conformer = ctx.conformer();

      const auto &rec_a = ctx.topology.atom(rec_idx_a);
      const auto &rec_b = ctx.topology.atom(rec_idx_b);

      if (are_residueids_close(mol, rec_a.atom, rec_b.atom, 1)) return InteractionType::None;
      if (!passes_hbond_distance_filter(mol, conformer, &rec_a.atom.get(), &rec_b.atom.get(), 0.1))         return InteractionType::None;
      if (!passes_hbond_angle_filter   (mol, conformer, &rec_a.atom.get(), &rec_b.atom.get(), 1.57, false)) return InteractionType::None; // 90 degrees
      return InteractionType::HydrogenBond;
    }
  };
}

ContactRecipe<AtomRec, AtomRec, WeakHbondParams> make_weak_hbond_recipe() {
  return {
    WeakHbondParams{},
    +[](const AtomRec &rec) { return (rec.type & AtomType::WeakHbondDonor) == AtomType::WeakHbondDonor; },
    +[](const AtomRec &rec) { return (rec.type & AtomType::HbondAcceptor)  == AtomType::HbondAcceptor; },
    +[](std::uint32_t rec_idx_a, std::uint32_t rec_idx_b, float dist_sq, const ContactContext& ctx) -> InteractionType {

      const auto &mol       = ctx.topology.molecule();
      const auto &conformer = ctx.conformer();

      const auto &rec_a = ctx.topology.atom(rec_idx_a);
      const auto &rec_b = ctx.topology.atom(rec_idx_b);

      if (are_residueids_close(mol, rec_a.atom, rec_b.atom, 1)) return InteractionType::None;
      if (!passes_hbond_distance_filter(mol, conformer, &rec_a.atom.get(), &rec_b.atom.get(), 0.1))        return InteractionType::None;
      if (!passes_hbond_angle_filter   (mol, conformer, &rec_a.atom.get(), &rec_b.atom.get(), 2.27, true)) return InteractionType::None; // 130 degrees
      return InteractionType::WeakHydrogenBond;
    }
  };
}

ContactRecipe<AtomRec, AtomRec, PolarHbondParams> make_polar_hbond_recipe() {
  return {
    PolarHbondParams{},
    +[](const AtomRec& rec) { return (rec.type & AtomType::HbondDonor)    == AtomType::HbondDonor; },
    +[](const AtomRec& rec) { return (rec.type & AtomType::HbondAcceptor) == AtomType::HbondAcceptor; },
    +[](std::uint32_t rec_idx_a, std::uint32_t rec_idx_b, float dist, const ContactContext& ctx) -> InteractionType {

      const auto &rec_a = ctx.topology.atom(rec_idx_a);
      const auto &rec_b = ctx.topology.atom(rec_idx_b);

      if (are_residueids_close(ctx.molecule(), rec_a.atom, rec_b.atom, 1)) return InteractionType::None;
      return InteractionType::PolarHydrogenBond;
    }
  };
}

ContactRecipe<AtomRec, AtomRec, WeakPolarHbondParams> make_weak_polar_hbond_recipe() {
  return {
    WeakPolarHbondParams{},
    +[](const AtomRec& rec) { return (rec.type & AtomType::WeakHbondDonor) == AtomType::WeakHbondDonor; },
    +[](const AtomRec& rec) { return (rec.type & AtomType::HbondAcceptor)  == AtomType::HbondAcceptor; },
    +[](std::uint32_t rec_idx_a, std::uint32_t rec_idx_b, float dist, const ContactContext& ctx) -> InteractionType {

      const auto& rec_a = ctx.topology.atom(rec_idx_a);
      const auto& rec_b = ctx.topology.atom(rec_idx_b);

      if (are_residueids_close(ctx.molecule(), rec_a.atom, rec_b.atom, 1)) return InteractionType::None;
      return InteractionType::WeakPolarHydrogenBond;
    }
  };
}

} // namespace lahuta::arpeggio
