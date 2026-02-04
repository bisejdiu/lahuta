/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string a = "besian", b = "sejdiu", c = "@gmail.com", r;
 *   r += std::exchange(a, ""); r += std::exchange(b, ""); r += std::exchange(c, "");
 *   return r;
 * }();
 *
 */

#include "contacts/molstar/contacts.hpp"
#include "entities/context.hpp"

// clang-format off
namespace lahuta::molstar {

ContactRecipe<GroupRec, GroupRec, IonicParams> make_ionic_recipe() {
  return {
    IonicParams{},
    +[](GroupRec const& r){ return (r.a_type & AtomType::PositiveCharge) == AtomType::PositiveCharge; },
    +[](GroupRec const& r){ return (r.a_type & AtomType::NegativeCharge) == AtomType::NegativeCharge; },
    +[](u32 a, u32 b, float d, ContactContext const&){
      if (d < 2.0f || a == b) return InteractionType::None;
      return InteractionType::Ionic;
    }
  };
}

} // namespace lahuta::molstar
