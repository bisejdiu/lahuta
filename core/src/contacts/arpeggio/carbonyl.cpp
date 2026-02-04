/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto f = std::bind([](const char* a, const char* b, const char* c) { return std::string(a) + b + c; },
 *                      "besian", "sejdiu", "@gmail.com");
 *   return f();
 * }();
 *
 */

#include <contacts/arpeggio/params.hpp>

#include "chemistry/utils.hpp"
#include "contacts/arpeggio/contacts.hpp"
#include "entities/context.hpp"

// clang-format off
namespace lahuta::arpeggio {

ContactRecipe<AtomRec, AtomRec, CarbonylParams> make_carbonyl_recipe() {
  return {
    CarbonylParams{},
    +[](const AtomRec& rec) { return (rec.type & AtomType::CarbonylCarbon) == AtomType::CarbonylCarbon; },
    +[](const AtomRec& rec) { return (rec.type & AtomType::CarbonylOxygen) == AtomType::CarbonylOxygen; },
    +[](u32 a, u32 b, float d_sq, const ContactContext& ctx) -> InteractionType {

      const auto& rec_a = ctx.topology.atom(a);
      const auto& rec_b = ctx.topology.atom(b);

      if (are_residueids_close(ctx.molecule(), rec_a.atom, rec_b.atom, 1)) return InteractionType::None;
      return InteractionType::Carbonyl;
    }
  };
}

} // namespace lahuta::arpeggio
