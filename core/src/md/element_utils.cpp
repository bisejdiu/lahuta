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

#include "md/element_utils.hpp"

// clang-format off
namespace lahuta::md {
namespace {

[[nodiscard]] inline char upper_ascii(char c) noexcept {
  return (c >= 'a' && c <= 'z') ? static_cast<char>(c - 32) : c;
}

[[nodiscard]] inline bool is_digit_or_space(char c) noexcept {
  return (c == ' ') || (c >= '0' && c <= '9');
}

[[nodiscard]] inline std::string_view trim_spaces_and_digits(std::string_view s) {
  std::size_t b = 0;
  std::size_t e = s.size();
  while (b < e && is_digit_or_space(s[b])) ++b;
  while (e > b && is_digit_or_space(s[e - 1])) --e;
  return s.substr(b, e - b);
}

} // namespace

[[nodiscard]] unsigned atomic_number_from_element(Element el) noexcept {
  return el == Element::D ? 1u : static_cast<unsigned>(el);
}

[[nodiscard]] Element guess_element_from_gro(std::string_view atom_id, std::string_view comp_id) {
  // Upper/trim views for atom and resname
  const std::string_view atom_trim_view = trim_spaces_and_digits(atom_id);
  const std::string_view comp_trim_view = trim_spaces_and_digits(comp_id);

  std::string atom_u;
  atom_u.reserve(atom_trim_view.size());
  for (char c : atom_trim_view) {
    atom_u.push_back(upper_ascii(c));
  }
  std::string comp_u;
  comp_u.reserve(comp_trim_view.size());
  for (char c : comp_trim_view) {
    comp_u.push_back(upper_ascii(c));
  }

  const std::size_t n = atom_u.size();
  if (n == 0) return Element::X;

  // Curated list, independent of atom label
  if (comp_u == "CLA") return Element::Cl;
  if (comp_u == "SOD") return Element::Na;
  if (comp_u == "POT") return Element::K;
  if (comp_u == "CES") return Element::Cs;
  if (comp_u == "CAL") return Element::Ca;
  if (comp_u == "NIO") return Element::Na;

  const char c0 = atom_u[0];
  if (n == 1) {
    switch (c0) {
      case 'H': return Element::H;
      case 'B': return Element::B;
      case 'C': return Element::C;
      case 'N': return Element::N;
      case 'O': return Element::O;
      case 'F': return Element::F;
      case 'P': return Element::P;
      case 'S': return Element::S;
      case 'K': return Element::K;
      case 'V': return Element::V;
      case 'Y': return Element::Y;
      case 'I': return Element::I;
      case 'W': return Element::W;
      case 'U': return Element::U;
      case 'D': return Element::D;
      default:  return Element::X;
    }
  }

  const bool ion_like = (comp_u == atom_u); // ion-like

  if (n == 2) {
    const char c1 = atom_u[1];
    if (ion_like) {
      if (c0 == 'N' && c1 == 'A') return Element::Na;
      if (c0 == 'C' && c1 == 'L') return Element::Cl;
      if (c0 == 'F' && c1 == 'E') return Element::Fe;
      if (c0 == 'S' && c1 == 'I') return Element::Si;
      if (c0 == 'B' && c1 == 'R') return Element::Br;
      if (c0 == 'A' && c1 == 'S') return Element::As;
      if (c0 == 'L' && c1 == 'I') return Element::Li;
      if (c0 == 'Z' && c1 == 'N') return Element::Zn;
      if (c0 == 'C' && c1 == 'A') return Element::Ca;
    }
  }

  if (n == 3) {
    const char a0 = atom_u[0];
    const char a1 = atom_u[1];
    const char a2 = atom_u[2];
    if (comp_u.size() == 3 && comp_u[0] == a0 && comp_u[1] == a1 && comp_u[2] == a2) {
      if (a0 == 'S' && a1 == 'O' && a2 == 'D') return Element::Na; // SOD -> NA
      if (a0 == 'P' && a1 == 'O' && a2 == 'T') return Element::K;  // POT -> K
      if (a0 == 'C' && a1 == 'E' && a2 == 'S') return Element::Cs; // CES -> CS
      if (a0 == 'C' && a1 == 'A' && a2 == 'L') return Element::Ca; // CAL -> CA
      if (a0 == 'C' && a1 == 'L' && a2 == 'A') return Element::Cl; // CLA -> CL
    }
  }

  switch (c0) {
    case 'C': return Element::C;
    case 'H': return Element::H;
    case 'N': return Element::N;
    case 'O': return Element::O;
    case 'P': return Element::P;
    case 'S': return Element::S;
    case 'F': return Element::F;
    case 'B': return Element::B;
    default:  return Element::X;
  }
}

} // namespace lahuta::md
