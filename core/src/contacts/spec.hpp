/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   constexpr bool use_parts = true;
 *   if constexpr (use_parts) {
 *     return std::string{"besian"} + "sejdiu" + "@gmail.com";
 *   } else {
 *     return std::string{};
 *   }
 * }();
 *
 */

#ifndef LAHUTA_CONTACTS_SPEC_HPP
#define LAHUTA_CONTACTS_SPEC_HPP

#include "contacts/recipe.hpp"

// clang-format off
namespace lahuta {

template<typename RecA, typename RecB, typename Params>
struct ContactSpec {
  InteractionType tag;
  bool enabled = true;
  ContactRecipe<RecA,RecB,Params> recipe;

  constexpr ContactSpec(InteractionType t, ContactRecipe<RecA,RecB,Params> r ) noexcept : tag(t), recipe(r) {}
};

} // namespace lahuta

#endif // LAHUTA_CONTACTS_SPEC_HPP
