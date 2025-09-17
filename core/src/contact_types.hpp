#pragma once

#include <string_view>

// clang-format off
namespace lahuta {

enum class AtomTypingMethod { Arpeggio, Molstar };

inline constexpr std::string_view contact_computer_name(AtomTypingMethod type) noexcept {
  switch (type) {
    case AtomTypingMethod::Arpeggio: return "Arpeggio";
    case AtomTypingMethod::Molstar:  return "MolStar";
  }
  return "unknown";
}

} // namespace lahuta
