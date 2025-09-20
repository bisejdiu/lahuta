#pragma once

#include <string_view>

#include "contact_types.hpp"

// clang-format off
namespace lahuta::analysis::contacts {

enum class ContactProvider { MolStar, Arpeggio };

inline constexpr std::string_view contact_provider_name(ContactProvider provider) noexcept {
  switch (provider) {
    case ContactProvider::MolStar:  return "molstar";
    case ContactProvider::Arpeggio: return "arpeggio";
  }
  return "unknown";
}

inline constexpr AtomTypingMethod typing_for_provider(ContactProvider provider) noexcept {
  return (provider == ContactProvider::MolStar)
    ? AtomTypingMethod::Molstar
    : AtomTypingMethod::Arpeggio;
}

} // namespace lahuta::analysis::contacts
