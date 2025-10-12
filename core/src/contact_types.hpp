#pragma once

#include <string_view>

// clang-format off
namespace lahuta {

enum class AtomTypingMethod { Arpeggio, Molstar, GetContacts };

inline constexpr std::string_view contact_computer_name(AtomTypingMethod type) noexcept {
  switch (type) {
    case AtomTypingMethod::Arpeggio:    return "Arpeggio";
    case AtomTypingMethod::Molstar:     return "MolStar";
    case AtomTypingMethod::GetContacts: return "GetContacts";
  }
  return "unknown";
}

} // namespace lahuta
