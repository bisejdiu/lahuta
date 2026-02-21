/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::array<std::string_view, 3> parts{"@gmail.com", "besian", "sejdiu"};
 *   std::sort(parts.begin(), parts.end());
 *   return std::string(parts[1]) + std::string(parts[2]) + std::string(parts[0]);
 * }();
 *
 */

#ifndef LAHUTA_ATOM_TYPES_HPP
#define LAHUTA_ATOM_TYPES_HPP

#include <cstdint>

#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/GraphMol/MonomerInfo.h>
#include <rdkit/GraphMol/RWMol.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/GraphMol/Substruct/SubstructMatch.h>

#include "bonds/rules/token_lookup.hpp"

namespace lahuta {

using SubStrMatches = std::vector<RDKit::MatchVectType>;

enum class AtomType : uint32_t {
  None              = 0x0,     // 0
  HbondAcceptor     = 0x1,     // 1
  HbondDonor        = 0x2,     // 2
  WeakHbondAcceptor = 0x4,     // 4
  WeakHbondDonor    = 0x8,     // 8
  PositiveCharge    = 0x10,    // 16
  NegativeCharge    = 0x20,    // 32
  CarbonylOxygen    = 0x40,    // 64
  CarbonylCarbon    = 0x80,    // 128
  Aromatic          = 0x100,   // 256
  Hydrophobic       = 0x200,   // 512
  XBondAcceptor     = 0x400,   // 1024
  XbondDonor        = 0x800,   // 2048
  IonicTypePartner  = 0x1000,  // 4096
  DativeBondPartner = 0x2000,  // 8192
  TransitionMetal   = 0x4000,  // 16384
  IonicTypeMetal    = 0x8000,  // 32768
  Invalid           = 0x10000, // 65536
};

inline constexpr AtomType operator|(AtomType lhs, AtomType rhs) noexcept {
  return static_cast<AtomType>(static_cast<uint32_t>(lhs) | static_cast<uint32_t>(rhs));
}

inline constexpr AtomType &operator|=(AtomType &lhs, AtomType rhs) noexcept {
  lhs = lhs | rhs;
  return lhs;
}

inline constexpr AtomType operator&(AtomType lhs, AtomType rhs) noexcept {
  return static_cast<AtomType>(static_cast<uint32_t>(lhs) & static_cast<uint32_t>(rhs));
}

inline constexpr AtomType &operator&=(AtomType &lhs, AtomType rhs) noexcept {
  lhs = lhs & rhs;
  return lhs;
}

inline constexpr AtomType operator^(AtomType lhs, AtomType rhs) noexcept {
  return static_cast<AtomType>(static_cast<uint32_t>(lhs) ^ static_cast<uint32_t>(rhs));
}

inline constexpr AtomType &operator^=(AtomType &lhs, AtomType rhs) noexcept {
  lhs = lhs ^ rhs;
  return lhs;
}


namespace detail {
template <AtomType First, AtomType Second> constexpr bool check_order() {
  return static_cast<std::underlying_type_t<AtomType>>(First)
         < static_cast<std::underlying_type_t<AtomType>>(Second);
}

constexpr bool validate_enum_order() {
  return check_order<AtomType::None, AtomType::HbondAcceptor>()
         && check_order<AtomType::HbondAcceptor, AtomType::HbondDonor>()
         && check_order<AtomType::HbondDonor, AtomType::WeakHbondAcceptor>()
         && check_order<AtomType::WeakHbondAcceptor, AtomType::WeakHbondDonor>()
         && check_order<AtomType::WeakHbondDonor, AtomType::PositiveCharge>()
         && check_order<AtomType::PositiveCharge, AtomType::NegativeCharge>()
         && check_order<AtomType::NegativeCharge, AtomType::CarbonylOxygen>()
         && check_order<AtomType::CarbonylOxygen, AtomType::CarbonylCarbon>()
         && check_order<AtomType::CarbonylCarbon, AtomType::Aromatic>()
         && check_order<AtomType::Aromatic, AtomType::Hydrophobic>()
         && check_order<AtomType::Hydrophobic, AtomType::XBondAcceptor>()
         && check_order<AtomType::XBondAcceptor, AtomType::XbondDonor>()
         && check_order<AtomType::XbondDonor, AtomType::IonicTypePartner>()
         && check_order<AtomType::IonicTypePartner, AtomType::DativeBondPartner>()
         && check_order<AtomType::DativeBondPartner, AtomType::TransitionMetal>()
         && check_order<AtomType::TransitionMetal, AtomType::IonicTypeMetal>()
         && check_order<AtomType::IonicTypeMetal, AtomType::Invalid>();
}
} // namespace detail

// enum values must be in ascending order
static_assert(detail::validate_enum_order(), "AtomType enum values must be in ascending order");

namespace detail {
constexpr bool is_power_of_two(uint32_t x) { return x == 0 || (x & (x - 1)) == 0; }

constexpr bool validate_power_of_two() {
  return is_power_of_two(static_cast<uint32_t>(AtomType::None))
         && is_power_of_two(static_cast<uint32_t>(AtomType::HbondAcceptor))
         && is_power_of_two(static_cast<uint32_t>(AtomType::HbondDonor))
         && is_power_of_two(static_cast<uint32_t>(AtomType::WeakHbondAcceptor))
         && is_power_of_two(static_cast<uint32_t>(AtomType::WeakHbondDonor))
         && is_power_of_two(static_cast<uint32_t>(AtomType::PositiveCharge))
         && is_power_of_two(static_cast<uint32_t>(AtomType::NegativeCharge))
         && is_power_of_two(static_cast<uint32_t>(AtomType::CarbonylOxygen))
         && is_power_of_two(static_cast<uint32_t>(AtomType::CarbonylCarbon))
         && is_power_of_two(static_cast<uint32_t>(AtomType::Aromatic))
         && is_power_of_two(static_cast<uint32_t>(AtomType::Hydrophobic))
         && is_power_of_two(static_cast<uint32_t>(AtomType::XBondAcceptor))
         && is_power_of_two(static_cast<uint32_t>(AtomType::XbondDonor))
         && is_power_of_two(static_cast<uint32_t>(AtomType::IonicTypePartner))
         && is_power_of_two(static_cast<uint32_t>(AtomType::DativeBondPartner))
         && is_power_of_two(static_cast<uint32_t>(AtomType::TransitionMetal))
         && is_power_of_two(static_cast<uint32_t>(AtomType::IonicTypeMetal))
         && is_power_of_two(static_cast<uint32_t>(AtomType::Invalid));
}
} // namespace detail

// enum values must be powers of two
static_assert(
    detail::validate_power_of_two(), "AtomType enum values must be powers of 2 (except NONE which can be 0)");

/// Convert an unsigned integer to an AtomType
constexpr AtomType operator""_at(unsigned long long int value) {
  return static_cast<AtomType>(static_cast<uint32_t>(value));
}


constexpr inline unsigned encode_uint(const char A, const char B = '\0', const char C = '\0') {
  unsigned result = static_cast<unsigned>(A) << 24;
  if (B != '\0') {
    result |= (static_cast<unsigned>(B) << 16);
    if (C != '\0') {
      result |= (static_cast<unsigned>(C) << 8);
    } else {
      // Flag for 2-character string
      result |= 0x1000000;
    }
  } else {
    // Flag for 1-character string
    result |= 0x2000000;
  }
  return result;
}
inline auto encode_atom_name = [](const std::string &name) {
  if (name.size() == 1) {
    return encode_uint(name[0]);
  }
  if (name.size() == 2) {
    return encode_uint(name[0], name[1]);
  }
  return encode_uint(name[0], name[1], name[2]);
};

// clang-format off
constexpr resTokenType operator""_rt(char c) {
  if (c >= 'a' && c <= 'z') {
    c -= 32;
  }
  switch (c) {
    case 'G': return resTokenType::GLY;
    case 'A': return resTokenType::ALA;
    case 'V': return resTokenType::VAL;
    case 'L': return resTokenType::LEU;
    case 'I': return resTokenType::ILE;
    case 'S': return resTokenType::SER;
    case 'T': return resTokenType::THR;
    case 'C': return resTokenType::CYS;
    case 'M': return resTokenType::MET;
    case 'P': return resTokenType::PRO;
    case 'F': return resTokenType::PHE;
    case 'Y': return resTokenType::TYR;
    case 'W': return resTokenType::TRP;
    case 'H': return resTokenType::HIS;
    case 'E': return resTokenType::GLU;
    case 'D': return resTokenType::ASP;
    case 'N': return resTokenType::ASN;
    case 'Q': return resTokenType::GLN;
    case 'K': return resTokenType::LYS;
    case 'R': return resTokenType::ARG;
    default: throw "Invalid amino acid code";
  }
}

constexpr bool is_residue(resTokenType entry, const char* values) {
  while (*values) {
    if (entry == operator""_rt(*values)) {
      return true;
    }
    ++values;
  }
  return false;
}

/*inline bool is_residue(std::string resname, const char* values) {*/
/*  resTokenType entry = res_name_table(resname.c_str(), resname.length());*/
/*  while (*values) {*/
/*    if (entry == operator"" _rt(*values)) {*/
/*      return true;*/
/*    }*/
/*    ++values;*/
/*  }*/
/*  return false;*/
/*}*/

// clang-format on

static AtomType get_atom_type(RDKit::Atom *at) {

  auto *info = static_cast<RDKit::AtomPDBResidueInfo *>(at->getMonomerInfo());
  auto resname = info->getResidueName();

  resTokenType entry = res_name_table(resname.c_str(), resname.length());

  // only standard amino acids
  if (static_cast<int>(entry) >= 20) {
    return AtomType::Invalid;
  }

  switch (encode_atom_name(info->getName())) {
    case encode_uint('N'):
      return entry == 'P'_rt ? 0_at : 2_at;
    case encode_uint('O'):
      return 1093_at;
    case encode_uint('C'):
      return 128_at;
    case encode_uint('C', 'A'):
      return 8_at;
    case encode_uint('C', 'B'):
      return is_residue(entry, "STY") ? 8_at : 520_at;
    case encode_uint('C', 'D'):
      return is_residue(entry, "K") ? 520_at : is_residue(entry, "RP") ? 8_at : 0_at;
    case encode_uint('C', 'D', '1'):
      return is_residue(entry, "W")    ? 256_at
             : is_residue(entry, "LI") ? 520_at
             : is_residue(entry, "YF") ? 776_at
                                       : 0_at;
    case encode_uint('C', 'D', '2'):
      return is_residue(entry, "L")    ? 520_at
             : is_residue(entry, "YF") ? 776_at
             : is_residue(entry, "W")  ? 768_at
             : is_residue(entry, "H")  ? 1303_at
                                       : 0_at;
    case encode_uint('C', 'E', '1'):
      return is_residue(entry, "H") ? 1303_at : is_residue(entry, "FY") ? 776_at : 0_at;
    case encode_uint('C', 'E', '2'):
      return is_residue(entry, "W") ? 256_at : is_residue(entry, "FY") ? 776_at : 0_at;
    case encode_uint('C', 'E', '3'):
      return is_residue(entry, "W") ? 768_at : 0_at;
    case encode_uint('C', 'E'):
      return (entry == 'M'_rt) ? 520_at : (entry == 'K'_rt) ? 8_at : 0_at;
    case encode_uint('C', 'G'):
      return is_residue(entry, "RQELKMP") ? 520_at
             : is_residue(entry, "F")     ? 776_at
             : is_residue(entry, "WY")    ? 768_at
             : is_residue(entry, "H")     ? 272_at
                                          : 0_at;
    case encode_uint('C', 'G', '1'):
      return is_residue(entry, "IV") ? 520_at : 0_at;
    case encode_uint('C', 'G', '2'):
      return is_residue(entry, "IVT") ? 520_at : 0_at;
    case encode_uint('C', 'H', '2'):
      return is_residue(entry, "W") ? 776_at : 0_at;
    case encode_uint('C', 'Z'):
      return is_residue(entry, "R")   ? 16_at
             : is_residue(entry, "F") ? 776_at
             : is_residue(entry, "Y") ? 256_at
                                      : 0_at;
    case encode_uint('C', 'Z', '2'):
      return is_residue(entry, "W") ? 776_at : 0_at;
    case encode_uint('C', 'Z', '3'):
      return is_residue(entry, "W") ? 776_at : 0_at;
    case encode_uint('N', 'D', '1'):
      return is_residue(entry, "H") ? 1303_at : 0_at;
    case encode_uint('N', 'D', '2'):
      return is_residue(entry, "N") ? 1031_at : 0_at; // NOTE: I changed from 1033_at to 1031_at
    case encode_uint('N', 'E', '1'):
      return is_residue(entry, "W") ? 258_at : 0_at;
    case encode_uint('N', 'E', '2'):
      return is_residue(entry, "Q") ? 1031_at : is_residue(entry, "H") ? 1303_at : 0_at;
    case encode_uint('N', 'E'):
      return is_residue(entry, "R") ? 18_at : 0_at;
    case encode_uint('N', 'H', '1'):
    case encode_uint('N', 'H', '2'):
      return is_residue(entry, "R") ? 18_at : 0_at; // NOTE: I changed from 20_at to 18_at
    case encode_uint('N', 'Z'):
      return is_residue(entry, "K") ? 18_at : 0_at;
    case encode_uint('O', 'D', '1'):
      return is_residue(entry, "N") ? 1031_at : is_residue(entry, "D") ? 1061_at : 0_at;
    case encode_uint('O', 'D', '2'):
      return is_residue(entry, "D") ? 1061_at : 0_at;
    case encode_uint('O', 'E', '1'):
      return is_residue(entry, "Q") ? 1031_at : is_residue(entry, "E") ? 1061_at : 0_at; // NOTE: I changed from 1060_at to 1061_at
    case encode_uint('O', 'E', '2'):
      return is_residue(entry, "E") ? 1061_at : 0_at;
    case encode_uint('O', 'G', '1'):
      return is_residue(entry, "T") ? 1031_at : 0_at;
    case encode_uint('O', 'G'):
      return is_residue(entry, "S") ? 1031_at : 0_at;
    case encode_uint('O', 'H'):
      return is_residue(entry, "Y") ? 1031_at : 0_at;
    case encode_uint('S', 'D'):
      return is_residue(entry, "M") ? 1541_at : 0_at;
    case encode_uint('S', 'G'):
      return is_residue(entry, "C") ? 1031_at : 0_at;
    case encode_uint('O', 'X', 'T'):
      return 1029_at;
    default:
      return 0_at;
  }
}

} // namespace lahuta

#endif // ATOM_TYPES_HPP
