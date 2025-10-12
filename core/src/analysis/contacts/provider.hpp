#pragma once

#include <string_view>

#include "contact_types.hpp"

// clang-format off
namespace lahuta::analysis::contacts {

enum class ContactProvider { MolStar, Arpeggio, GetContacts };

inline constexpr std::string_view contact_provider_name(ContactProvider provider) noexcept {
  switch (provider) {
    case ContactProvider::MolStar:     return "molstar";
    case ContactProvider::Arpeggio:    return "arpeggio";
    case ContactProvider::GetContacts: return "getcontacts";
  }
  return "unknown";
}

inline constexpr AtomTypingMethod typing_for_provider(ContactProvider provider) noexcept {
  switch (provider) {
    case ContactProvider::Arpeggio:    return AtomTypingMethod::Arpeggio;
    case ContactProvider::GetContacts: return AtomTypingMethod::GetContacts;
    case ContactProvider::MolStar:
    default:
      return AtomTypingMethod::Molstar;
  }
}

} // namespace lahuta::analysis::contacts
