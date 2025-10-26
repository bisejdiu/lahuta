#ifndef LAHUTA_MODELS_DSSP_HPP
#define LAHUTA_MODELS_DSSP_HPP

#include <cstdint>

// clang-format off
namespace lahuta {

//
// This enum currently reflects AlphaFold mmCIF `_struct_conf` annotations
// rather than canonical DSSP output. I plan to generalize this to true
// DSSP semantics (H, B, E, G, I, T, S, Coil).  - Besian, October 2025
//

enum class DSSPAssignment : std::uint8_t {
  Coil = 0,
  AlphaHelix,        // H
  Helix3_10,         // G
  HelixPi,           // I
  PolyProlineHelix,  // P (left-handed polyproline)
  Strand,            // E
  Turn,              // T
  Bend               // S
};

static_assert(sizeof(DSSPAssignment) == sizeof(std::uint8_t), "DSSPAssignment must remain 1 byte");

} // namespace lahuta

#endif // LAHUTA_MODELS_DSSP_HPP
