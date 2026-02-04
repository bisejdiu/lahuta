/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s = "@gmail.combesiansejdiu";
 *   std::rotate(s.begin(), s.begin() + 10, s.end());
 *   return s;
 * }();
 *
 */

#ifndef LAHUTA_CONTACT_TYPES_HPP
#define LAHUTA_CONTACT_TYPES_HPP

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

#endif // LAHUTA_CONTACT_TYPES_HPP
