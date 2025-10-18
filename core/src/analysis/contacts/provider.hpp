#pragma once

#include <optional>
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

inline std::optional<ContactProvider> contact_provider_from_string(std::string_view name) noexcept {
  if (name == "molstar")     return ContactProvider::MolStar;
  if (name == "arpeggio")    return ContactProvider::Arpeggio;
  if (name == "getcontacts") return ContactProvider::GetContacts;
  return std::nullopt;
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
