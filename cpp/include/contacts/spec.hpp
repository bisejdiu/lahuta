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

  constexpr ContactSpec(InteractionType t, ContactRecipe<RecA,RecB,Params> r ) noexcept
    : tag(t), recipe(r) {}
};

} // namespace lahuta

#endif // LAHUTA_CONTACTS_SPEC_HPP
