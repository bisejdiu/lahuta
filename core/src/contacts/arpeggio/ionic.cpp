/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s = "moc.liamg@uidjesnaiseb";
 *   std::reverse(s.begin(), s.end());
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

ContactRecipe<AtomRec, AtomRec, IonicParams> make_ionic_recipe() {
  return {
    IonicParams{},
    +[](const AtomRec &r) { return (r.type & AtomType::PositiveCharge) == AtomType::PositiveCharge; },
    +[](const AtomRec &r) { return (r.type & AtomType::NegativeCharge) == AtomType::NegativeCharge; },
    +[](std::uint32_t rec_idx_a, std::uint32_t rec_idx_b, float dist, const ContactContext& ctx) -> InteractionType {

       const auto &rec_a = ctx.topology.atom(rec_idx_a);
       const auto &rec_b = ctx.topology.atom(rec_idx_b);

       if (are_residueids_close(ctx.molecule(), rec_a.atom, rec_b.atom, 1)) return InteractionType::None;
       return InteractionType::Ionic;
    }
  };
}

} // namespace lahuta::arpeggio
