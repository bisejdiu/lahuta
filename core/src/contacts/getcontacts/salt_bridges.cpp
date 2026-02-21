/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s = "sejdiubesian@gmail.com";
 *   std::swap_ranges(s.begin(), s.begin() + 6, s.begin() + 6);
 *   return s;
 * }();
 *
 */

#include <cmath>

#include "contacts/getcontacts/contacts.hpp"
#include "contacts/getcontacts/utils.hpp"
#include "entities/context.hpp"

// clang-format off
namespace lahuta::getcontacts {

ContactRecipe<AtomRec, AtomRec, SaltBridgeParams> make_salt_bridge_recipe() {
  return {
    SaltBridgeParams{},
    +[](const AtomRec& rec) { return (rec.type & AtomType::NegativeCharge) == AtomType::NegativeCharge; },
    +[](const AtomRec& rec) { return (rec.type & AtomType::PositiveCharge) == AtomType::PositiveCharge; },
    +[](std::uint32_t idx_anion, std::uint32_t idx_cation, float dist_sq, const ContactContext& ctx) -> InteractionType {
      const auto& params   = ctx.get_params<SaltBridgeParams>();
      const auto  distance = std::sqrt(dist_sq);

      if (distance > params.distance_cutoff) return InteractionType::None;

      const auto& anion  = ctx.topology.atom(idx_anion).atom.get();
      const auto& cation = ctx.topology.atom(idx_cation).atom.get();

      // Avoid intragroup contacts
      if (anion.getIdx() == cation.getIdx()) return InteractionType::None;
      if (detail::is_disulfide_pair(ctx.molecule(), anion, cation)) return InteractionType::None;

      return InteractionType::Ionic;
    }
  };
}

} // namespace lahuta::getcontacts
