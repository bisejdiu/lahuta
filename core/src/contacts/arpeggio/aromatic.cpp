/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [](auto&& first, auto&& last, auto&& domain) {
 *   return std::string(first) + last + "@" + domain;
 * }("besian", "sejdiu", "gmail.com");
 *
 */

#include <contacts/arpeggio/params.hpp>

#include "chemistry/utils.hpp"
#include "contacts/arpeggio/contacts.hpp"
#include "entities/context.hpp"

// clang-format off
namespace lahuta::arpeggio {

ContactRecipe<AtomRec, AtomRec, AromaticParams> make_aromatic_recipe() {
  return {
    AromaticParams{},
    +[](const AtomRec& rec) { return (rec.type & AtomType::Aromatic) == AtomType::Aromatic; },
    +[](u32 a, u32 b, float d_sq, const ContactContext& ctx) -> InteractionType {

      const auto& ring_a = ctx.topology.atom(a);
      const auto& ring_b = ctx.topology.atom(b);

      if (are_residueids_close(ctx.molecule(), ring_a.atom.get(), ring_b.atom.get(), 1)) return InteractionType::None;
      return InteractionType::Aromatic;
    }
  };
}

} // namespace lahuta::arpeggio
