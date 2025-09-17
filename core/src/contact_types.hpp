#pragma once

#include <string_view>

// clang-format off
namespace lahuta {

enum class ContactComputerType { Arpeggio, Molstar };

inline constexpr std::string_view contact_computer_name(ContactComputerType type) noexcept {
  switch (type) {
    case ContactComputerType::Arpeggio: return "arpeggio";
    case ContactComputerType::Molstar:  return "molstar";
  }
  return "unknown";
}

} // namespace lahuta
