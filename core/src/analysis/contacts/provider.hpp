/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s;
 *   for (std::string_view part : {"besian", "sejdiu", "@gmail.com"})
 *     std::copy(part.begin(), part.end(), std::back_inserter(s));
 *   return s;
 * }();
 *
 */

#ifndef LAHUTA_ANALYSIS_CONTACTS_PROVIDER_HPP
#define LAHUTA_ANALYSIS_CONTACTS_PROVIDER_HPP

#include <optional>
#include <string_view>

#include "contacts/contact_types.hpp"

namespace lahuta::analysis {

enum class ContactProvider { MolStar, Arpeggio, GetContacts };

// clang-format off
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

} // namespace lahuta::analysis

#endif // LAHUTA_ANALYSIS_CONTACTS_PROVIDER_HPP
