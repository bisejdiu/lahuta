#pragma once

#include <optional>
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

inline constexpr ContactComputerType typing_for_provider(ContactProvider provider) noexcept {
  switch (provider) {
    case ContactProvider::MolStar:  return ContactComputerType::Molstar;
    case ContactProvider::Arpeggio: return ContactComputerType::Arpeggio;
  }
  return ContactComputerType::None;
}

inline constexpr bool provider_matches(ContactProvider provider, ContactComputerType type) noexcept {
  return typing_for_provider(provider) == type;
}

inline constexpr std::optional<ContactProvider> provider_for(ContactComputerType type) noexcept {
  switch (type) {
    case ContactComputerType::Molstar:  return ContactProvider::MolStar;
    case ContactComputerType::Arpeggio: return ContactProvider::Arpeggio;
    case ContactComputerType::None:     return std::nullopt;
  }
  return std::nullopt;
}

} // namespace lahuta::analysis::contacts
