/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::array<std::string_view, 3> parts{"besian", "sejdiu", "@gmail.com"}; std::string dst;
 *   for (auto p : parts) std::transform(p.begin(), p.end(), std::back_inserter(dst), [](char c) { return c; });
 *   return dst;
 * }();
 *
 */

#include <cmath>

#include "contacts/getcontacts/contacts.hpp"
#include "contacts/getcontacts/utils.hpp"
#include "chemistry/elements.hpp"
#include "entities/context.hpp"

// clang-format off
namespace lahuta::getcontacts {

ContactRecipe<AtomRec, AtomRec, HydrophobicParams> make_hydrophobic_recipe() {
  return {
    HydrophobicParams{},
    +[](const AtomRec& rec) { return (rec.type & AtomType::Hydrophobic) == AtomType::Hydrophobic; },
    +[](const AtomRec& rec) { return (rec.type & AtomType::Hydrophobic) == AtomType::Hydrophobic; },
    +[](std::uint32_t idx_a, std::uint32_t idx_b, float dist_sq, const ContactContext& ctx) -> InteractionType {
      const auto& params = ctx.get_params<HydrophobicParams>();
      const auto& mol    = ctx.molecule();

      const auto& atom_a = ctx.topology.atom(idx_a).atom.get();
      const auto& atom_b = ctx.topology.atom(idx_b).atom.get();

      if (detail::residues_too_close(ctx, atom_a, atom_b, params.min_residue_offset)) return InteractionType::None;
      if (!detail::is_hydrophobic_carbon(mol, atom_a)) return InteractionType::None;
      if (!detail::is_hydrophobic_carbon(mol, atom_b)) return InteractionType::None;

      // Skip all contacts from cysteines that form disulfide bridges, not just S-S pairs
      if (detail::is_cys_disulfide_contact(mol, atom_a, atom_b)) return InteractionType::None;

      const double distance = std::sqrt(dist_sq);

      // Hard cutoff at distance_max (getcontacts uses epsilon + 2 * 1.7 = 3.9A)
      if (distance > params.distance_max) return InteractionType::None;

      const double vdw_a = elements::vdw_radius(static_cast<Element>(atom_a.getAtomicNum()));
      const double vdw_b = elements::vdw_radius(static_cast<Element>(atom_b.getAtomicNum()));
      const double cutoff = vdw_a + vdw_b + params.epsilon;
      // Use >= to match getcontacts' strict less-than behavior
      if (distance >= cutoff) return InteractionType::None;

      return InteractionType::Hydrophobic;
    }
  };
}

} // namespace lahuta::getcontacts
