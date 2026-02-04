/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   constexpr auto f = []() constexpr { return "besian"; };
 *   constexpr auto l = []() constexpr { return "sejdiu"; };
 *   constexpr auto d = []() constexpr { return "@gmail.com"; };
 *   return std::string(f()) + l() + d();
 * }();
 *
 */

#ifndef LAHUTA_ANALYSIS_SASA_PROTOR_RADII_HPP
#define LAHUTA_ANALYSIS_SASA_PROTOR_RADII_HPP

#include <array>
#include <cstddef>
#include <cstdint>
#include <string_view>

namespace lahuta::analysis {

// ProtOr radii (Tsai et al. 1999), matches FreeSASA/RustSASA defaults.
struct ProtorAtomRadius {
  const char *name;
  double radius;
};

// clang-format off
inline constexpr std::array<ProtorAtomRadius, 4> ProtorGly = {{
    {"N",   1.64},
    {"CA",  1.88},
    {"C",   1.61},
    {"O",   1.42},
}};

inline constexpr std::array<ProtorAtomRadius, 5> ProtorAla = {{
    {"N",   1.64},
    {"CA",  1.88},
    {"C",   1.61},
    {"CB",  1.88},
    {"O",   1.42},
}};

inline constexpr std::array<ProtorAtomRadius, 7> ProtorVal = {{
    {"N",   1.64},
    {"CA",  1.88},
    {"C",   1.61},
    {"CB",  1.88},
    {"O",   1.42},
    {"CG1", 1.88},
    {"CG2", 1.88},
}};

inline constexpr std::array<ProtorAtomRadius, 8> ProtorLeu = {{
    {"N",   1.64},
    {"CA",  1.88},
    {"C",   1.61},
    {"CB",  1.88},
    {"O",   1.42},
    {"CG",  1.88},
    {"CD1", 1.88},
    {"CD2", 1.88},
}};

inline constexpr std::array<ProtorAtomRadius, 8> ProtorIle = {{
    {"N",   1.64},
    {"CA",  1.88},
    {"C",   1.61},
    {"CB",  1.88},
    {"O",   1.42},
    {"CG1", 1.88},
    {"CG2", 1.88},
    {"CD1", 1.88},
}};

inline constexpr std::array<ProtorAtomRadius, 6> ProtorSer = {{
    {"N",   1.64},
    {"CA",  1.88},
    {"C",   1.61},
    {"CB",  1.88},
    {"O",   1.42},
    {"OG",  1.46},
}};

inline constexpr std::array<ProtorAtomRadius, 7> ProtorThr = {{
    {"N",   1.64},
    {"CA",  1.88},
    {"C",   1.61},
    {"CB",  1.88},
    {"O",   1.42},
    {"CG2", 1.88},
    {"OG1", 1.46},
}};

inline constexpr std::array<ProtorAtomRadius, 6> ProtorCys = {{
    {"N",   1.64},
    {"CA",  1.88},
    {"C",   1.61},
    {"CB",  1.88},
    {"O",   1.42},
    {"SG",  1.77},
}};

inline constexpr std::array<ProtorAtomRadius, 8> ProtorMet = {{
    {"N",   1.64},
    {"CA",  1.88},
    {"C",   1.61},
    {"CB",  1.88},
    {"O",   1.42},
    {"CG",  1.88},
    {"SD",  1.77},
    {"CE",  1.88},
}};

inline constexpr std::array<ProtorAtomRadius, 7> ProtorPro = {{
    {"N",   1.64},
    {"CA",  1.88},
    {"C",   1.61},
    {"CB",  1.88},
    {"O",   1.42},
    {"CG",  1.88},
    {"CD",  1.88},
}};

inline constexpr std::array<ProtorAtomRadius, 11> ProtorPhe = {{
    {"N",   1.64},
    {"CA",  1.88},
    {"C",   1.61},
    {"CB",  1.88},
    {"O",   1.42},
    {"CG",  1.61},
    {"CD1", 1.76},
    {"CD2", 1.76},
    {"CE1", 1.76},
    {"CE2", 1.76},
    {"CZ",  1.76},
}};

inline constexpr std::array<ProtorAtomRadius, 12> ProtorTyr = {{
    {"N",   1.64},
    {"CA",  1.88},
    {"C",   1.61},
    {"CB",  1.88},
    {"O",   1.42},
    {"CG",  1.61},
    {"CD1", 1.76},
    {"CD2", 1.76},
    {"CE1", 1.76},
    {"CE2", 1.76},
    {"OH",  1.46},
    {"CZ",  1.61},
}};

inline constexpr std::array<ProtorAtomRadius, 14> ProtorTrp = {{
    {"N",   1.64},
    {"CA",  1.88},
    {"C",   1.61},
    {"CB",  1.88},
    {"O",   1.42},
    {"CG",  1.61},
    {"CD1", 1.76},
    {"CD2", 1.61},
    {"CE2", 1.61},
    {"CE3", 1.76},
    {"NE1", 1.64},
    {"CH2", 1.76},
    {"CZ2", 1.76},
    {"CZ3", 1.76},
}};

inline constexpr std::array<ProtorAtomRadius, 10> ProtorHis = {{
    {"N",   1.64},
    {"CA",  1.88},
    {"C",   1.61},
    {"CB",  1.88},
    {"O",   1.42},
    {"CG",  1.61},
    {"CD2", 1.76},
    {"ND1", 1.64},
    {"CE1", 1.76},
    {"NE2", 1.64},
}};

inline constexpr std::array<ProtorAtomRadius, 9> ProtorGlu = {{
    {"N",   1.64},
    {"CA",  1.88},
    {"C",   1.61},
    {"CB",  1.88},
    {"O",   1.42},
    {"CG",  1.88},
    {"CD",  1.61},
    {"OE1", 1.42},
    {"OE2", 1.46},
}};

inline constexpr std::array<ProtorAtomRadius, 8> ProtorAsp = {{
    {"N",   1.64},
    {"CA",  1.88},
    {"C",   1.61},
    {"CB",  1.88},
    {"O",   1.42},
    {"CG",  1.61},
    {"OD1", 1.42},
    {"OD2", 1.46},
}};

inline constexpr std::array<ProtorAtomRadius, 8> ProtorAsn = {{
    {"N",   1.64},
    {"CA",  1.88},
    {"C",   1.61},
    {"CB",  1.88},
    {"O",   1.42},
    {"CG",  1.61},
    {"ND2", 1.64},
    {"OD1", 1.42},
}};

inline constexpr std::array<ProtorAtomRadius, 9> ProtorGln = {{
    {"N",   1.64},
    {"CA",  1.88},
    {"C",   1.61},
    {"CB",  1.88},
    {"O",   1.42},
    {"CG",  1.88},
    {"CD",  1.61},
    {"NE2", 1.64},
    {"OE1", 1.42},
}};

inline constexpr std::array<ProtorAtomRadius, 9> ProtorLys = {{
    {"N",   1.64},
    {"CA",  1.88},
    {"C",   1.61},
    {"CB",  1.88},
    {"O",   1.42},
    {"CG",  1.88},
    {"CD",  1.88},
    {"CE",  1.88},
    {"NZ",  1.64},
}};

inline constexpr std::array<ProtorAtomRadius, 11> ProtorArg = {{
    {"N",   1.64},
    {"CA",  1.88},
    {"C",   1.61},
    {"CB",  1.88},
    {"O",   1.42},
    {"CG",  1.88},
    {"CD",  1.88},
    {"NE",  1.64},
    {"NH1", 1.64},
    {"NH2", 1.64},
    {"CZ",  1.61},
}};

inline constexpr double ProtorOxtRadius = 1.46;

[[nodiscard]] constexpr std::uint32_t make_residue_key(char a, char b, char c) noexcept {
  return (static_cast<std::uint32_t>(static_cast<unsigned char>(a)) << 16) |
         (static_cast<std::uint32_t>(static_cast<unsigned char>(b)) << 8) |
          static_cast<std::uint32_t>(static_cast<unsigned char>(c));
}

// Look up ProtOr radius for a given residue and atom name.
// Returns -1.0 if not found.
[[nodiscard]] inline double protor_radius_for_atom(std::string_view residue_name, std::string_view atom_name) {
  const ProtorAtomRadius *atoms = nullptr;
  std::size_t count             = 0;

  if (residue_name.size() != 3) return -1.0;

  switch (make_residue_key(residue_name[0], residue_name[1], residue_name[2])) {
    case make_residue_key('G', 'L', 'Y'):
      atoms = ProtorGly.data();
      count = ProtorGly.size();
      break;
    case make_residue_key('A', 'L', 'A'):
      atoms = ProtorAla.data();
      count = ProtorAla.size();
      break;
    case make_residue_key('V', 'A', 'L'):
      atoms = ProtorVal.data();
      count = ProtorVal.size();
      break;
    case make_residue_key('L', 'E', 'U'):
      atoms = ProtorLeu.data();
      count = ProtorLeu.size();
      break;
    case make_residue_key('I', 'L', 'E'):
      atoms = ProtorIle.data();
      count = ProtorIle.size();
      break;
    case make_residue_key('S', 'E', 'R'):
      atoms = ProtorSer.data();
      count = ProtorSer.size();
      break;
    case make_residue_key('T', 'H', 'R'):
      atoms = ProtorThr.data();
      count = ProtorThr.size();
      break;
    case make_residue_key('C', 'Y', 'S'):
      atoms = ProtorCys.data();
      count = ProtorCys.size();
      break;
    case make_residue_key('M', 'E', 'T'):
      atoms = ProtorMet.data();
      count = ProtorMet.size();
      break;
    case make_residue_key('P', 'R', 'O'):
      atoms = ProtorPro.data();
      count = ProtorPro.size();
      break;
    case make_residue_key('P', 'H', 'E'):
      atoms = ProtorPhe.data();
      count = ProtorPhe.size();
      break;
    case make_residue_key('T', 'Y', 'R'):
      atoms = ProtorTyr.data();
      count = ProtorTyr.size();
      break;
    case make_residue_key('T', 'R', 'P'):
      atoms = ProtorTrp.data();
      count = ProtorTrp.size();
      break;
    case make_residue_key('H', 'I', 'S'):
      atoms = ProtorHis.data();
      count = ProtorHis.size();
      break;
    case make_residue_key('G', 'L', 'U'):
      atoms = ProtorGlu.data();
      count = ProtorGlu.size();
      break;
    case make_residue_key('A', 'S', 'P'):
      atoms = ProtorAsp.data();
      count = ProtorAsp.size();
      break;
    case make_residue_key('A', 'S', 'N'):
      atoms = ProtorAsn.data();
      count = ProtorAsn.size();
      break;
    case make_residue_key('G', 'L', 'N'):
      atoms = ProtorGln.data();
      count = ProtorGln.size();
      break;
    case make_residue_key('L', 'Y', 'S'):
      atoms = ProtorLys.data();
      count = ProtorLys.size();
      break;
    case make_residue_key('A', 'R', 'G'):
      atoms = ProtorArg.data();
      count = ProtorArg.size();
      break;
    default:
      return -1.0;
  }

  for (std::size_t i = 0; i < count; ++i) {
    if (atom_name == atoms[i].name) {
      return atoms[i].radius;
    }
  }
  return -1.0;
}

} // namespace lahuta::analysis

#endif // LAHUTA_ANALYSIS_SASA_PROTOR_RADII_HPP
