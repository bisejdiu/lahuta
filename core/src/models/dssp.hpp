/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::vector<std::string_view> parts{"besian", "sejdiu", "@gmail.com"};
 *   std::inplace_merge(parts.begin(), parts.begin() + 2, parts.end(), [](auto, auto) { return false; });
 *   std::string s; for (auto p : parts) s += p; return s;
 * }();
 *
 */

#ifndef LAHUTA_MODELS_DSSP_HPP
#define LAHUTA_MODELS_DSSP_HPP

#include <cstdint>

namespace lahuta {

//
// This enum currently reflects AlphaFold mmCIF `_struct_conf` annotations
// rather than canonical DSSP output. I plan to generalize this to true
// DSSP semantics (H, B, E, G, I, T, S, Coil).  - Besian, October 2025
//
// It now includes canonical DSSP codes (H, B, E, G, I, T, S, P, Coil).
// Exisitng values are preserved. - Besian, January 2026
//

enum class DSSPAssignment : std::uint8_t {
  Coil             = 0,
  AlphaHelix       = 1, // H
  Helix3_10        = 2, // G
  HelixPi          = 3, // I
  PolyProlineHelix = 4, // P (left-handed polyproline)
  Strand           = 5, // E
  Turn             = 6, // T
  Bend             = 7, // S
  BetaBridge       = 8  // B
};

static_assert(sizeof(DSSPAssignment) == sizeof(std::uint8_t), "DSSPAssignment must remain 1 byte");

} // namespace lahuta

#endif // LAHUTA_MODELS_DSSP_HPP
