#ifndef LAHUTA_ATOM_TYPES_HPP
#define LAHUTA_ATOM_TYPES_HPP

#include <cstdint>
#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/GraphMol/MonomerInfo.h>
#include <rdkit/GraphMol/RWMol.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/GraphMol/Substruct/SubstructMatch.h>

#include "bonds/token-gperf-generated.hpp"
#include "bonds/token.h"

namespace lahuta {

using SubStrMatches = std::vector<RDKit::MatchVectType>;

enum class AtomType : uint32_t {
  NONE = 0x0,                 // 0
  HBOND_ACCEPTOR = 0x1,       // 1
  HBOND_DONOR = 0x2,          // 2
  WEAK_HBOND_ACCEPTOR = 0x4,  // 4
  WEAK_HBOND_DONOR = 0x8,     // 8
  POS_IONISABLE = 0x10,       // 16
  NEG_IONISABLE = 0x20,       // 32
  CARBONYL_OXYGEN = 0x40,     // 64
  CARBONYL_CARBON = 0x80,     // 128
  AROMATIC = 0x100,           // 256
  HYDROPHOBIC = 0x200,        // 512
  XBOND_ACCEPTOR = 0x400,     // 1024
  XBOND_DONOR = 0x800,        // 2048
  IonicTypePartner = 0x1000,  // 4096
  DativeBondPartner = 0x2000, // 8192
  TransitionMetal = 0x4000,   // 16384
  IonicTypeMetal = 0x8000,    // 32768
  INVALID = 0x10000,          // 65536
};

namespace detail {
template <AtomType First, AtomType Second> constexpr bool check_order() {
  return static_cast<std::underlying_type_t<AtomType>>(First)
         < static_cast<std::underlying_type_t<AtomType>>(Second);
}

constexpr bool validate_enum_order() {
  // Check that each enum value is less than the next one
  return check_order<AtomType::NONE, AtomType::HBOND_ACCEPTOR>()
         && check_order<AtomType::HBOND_ACCEPTOR, AtomType::HBOND_DONOR>()
         && check_order<AtomType::HBOND_DONOR, AtomType::WEAK_HBOND_ACCEPTOR>()
         && check_order<AtomType::WEAK_HBOND_ACCEPTOR, AtomType::WEAK_HBOND_DONOR>()
         && check_order<AtomType::WEAK_HBOND_DONOR, AtomType::POS_IONISABLE>()
         && check_order<AtomType::POS_IONISABLE, AtomType::NEG_IONISABLE>()
         && check_order<AtomType::NEG_IONISABLE, AtomType::CARBONYL_OXYGEN>()
         && check_order<AtomType::CARBONYL_OXYGEN, AtomType::CARBONYL_CARBON>()
         && check_order<AtomType::CARBONYL_CARBON, AtomType::AROMATIC>()
         && check_order<AtomType::AROMATIC, AtomType::HYDROPHOBIC>()
         && check_order<AtomType::HYDROPHOBIC, AtomType::XBOND_ACCEPTOR>()
         && check_order<AtomType::XBOND_ACCEPTOR, AtomType::XBOND_DONOR>()
         && check_order<AtomType::XBOND_DONOR, AtomType::IonicTypePartner>()
         && check_order<AtomType::IonicTypePartner, AtomType::DativeBondPartner>()
         && check_order<AtomType::DativeBondPartner, AtomType::TransitionMetal>()
         && check_order<AtomType::TransitionMetal, AtomType::IonicTypeMetal>()
         && check_order<AtomType::IonicTypeMetal, AtomType::INVALID>();
}
} // namespace detail

// enum values must be in ascending order
static_assert(detail::validate_enum_order(), "AtomType enum values must be in ascending order");

namespace detail {
constexpr bool is_power_of_two(uint32_t x) { return x == 0 || (x & (x - 1)) == 0; }

constexpr bool validate_power_of_two() {
  return is_power_of_two(static_cast<uint32_t>(AtomType::NONE))
         && is_power_of_two(static_cast<uint32_t>(AtomType::HBOND_ACCEPTOR))
         && is_power_of_two(static_cast<uint32_t>(AtomType::HBOND_DONOR))
         && is_power_of_two(static_cast<uint32_t>(AtomType::WEAK_HBOND_ACCEPTOR))
         && is_power_of_two(static_cast<uint32_t>(AtomType::WEAK_HBOND_DONOR))
         && is_power_of_two(static_cast<uint32_t>(AtomType::POS_IONISABLE))
         && is_power_of_two(static_cast<uint32_t>(AtomType::NEG_IONISABLE))
         && is_power_of_two(static_cast<uint32_t>(AtomType::CARBONYL_OXYGEN))
         && is_power_of_two(static_cast<uint32_t>(AtomType::CARBONYL_CARBON))
         && is_power_of_two(static_cast<uint32_t>(AtomType::AROMATIC))
         && is_power_of_two(static_cast<uint32_t>(AtomType::HYDROPHOBIC))
         && is_power_of_two(static_cast<uint32_t>(AtomType::XBOND_ACCEPTOR))
         && is_power_of_two(static_cast<uint32_t>(AtomType::XBOND_DONOR))
         && is_power_of_two(static_cast<uint32_t>(AtomType::IonicTypePartner))
         && is_power_of_two(static_cast<uint32_t>(AtomType::DativeBondPartner))
         && is_power_of_two(static_cast<uint32_t>(AtomType::TransitionMetal))
         && is_power_of_two(static_cast<uint32_t>(AtomType::IonicTypeMetal))
         && is_power_of_two(static_cast<uint32_t>(AtomType::INVALID));
}
} // namespace detail

// enum values must be powers of two
static_assert(
    detail::validate_power_of_two(), "AtomType enum values must be powers of 2 (except NONE which can be 0)");

/// Convert an unsigned integer to an AtomType
constexpr AtomType operator""_at(unsigned long long int value) {
  return static_cast<AtomType>(static_cast<uint32_t>(value));
}

inline AtomType string_to_atom_type(const std::string &flag_name) {
  static const std::unordered_map<std::string, AtomType> stringToEnum = {
      {"NONE", AtomType::NONE},
      {"HBOND_ACCEPTOR", AtomType::HBOND_ACCEPTOR},
      {"HBOND_DONOR", AtomType::HBOND_DONOR},
      {"WEAK_HBOND_ACCEPTOR", AtomType::WEAK_HBOND_ACCEPTOR},
      {"WEAK_HBOND_DONOR", AtomType::WEAK_HBOND_DONOR},
      {"POS_IONISABLE", AtomType::POS_IONISABLE},
      {"NEG_IONISABLE", AtomType::NEG_IONISABLE},
      {"CARBONYL_OXYGEN", AtomType::CARBONYL_OXYGEN},
      {"CARBONYL_CARBON", AtomType::CARBONYL_CARBON},
      {"AROMATIC", AtomType::AROMATIC},
      {"HYDROPHOBIC", AtomType::HYDROPHOBIC},
      {"XBOND_ACCEPTOR", AtomType::XBOND_ACCEPTOR},
      {"XBOND_DONOR", AtomType::XBOND_DONOR},
      {"IONIC_TYPE_PARTNER", AtomType::IonicTypePartner},
      {"DATIVE_PARTNER", AtomType::DativeBondPartner},
      {"TRANSITION_METAL", AtomType::TransitionMetal},
      {"IONIC_TYPE_METAL", AtomType::IonicTypeMetal},
      {"INVALID", AtomType::INVALID}};

  auto it = stringToEnum.find(flag_name);
  if (it != stringToEnum.end()) {
    return it->second;
  }
  throw std::invalid_argument("Invalid AtomType flag name: " + flag_name);
}

constexpr AtomType operator|(AtomType lhs, AtomType rhs) {
  return static_cast<AtomType>(static_cast<uint32_t>(lhs) | static_cast<uint32_t>(rhs));
}

constexpr AtomType &operator|=(AtomType &lhs, AtomType rhs) {
  lhs = lhs | rhs;
  return lhs;
}

constexpr AtomType operator&(AtomType lhs, AtomType rhs) {
  return static_cast<AtomType>(static_cast<uint32_t>(lhs) & static_cast<uint32_t>(rhs));
}

constexpr AtomType &operator&=(AtomType &lhs, AtomType rhs) {
  lhs = lhs & rhs;
  return lhs;
}

constexpr AtomType operator^(AtomType lhs, AtomType rhs) {
  return static_cast<AtomType>(static_cast<uint32_t>(lhs) ^ static_cast<uint32_t>(rhs));
}

constexpr AtomType &operator^=(AtomType &lhs, AtomType rhs) {
  lhs = lhs ^ rhs;
  return lhs;
}

constexpr AtomType operator~(AtomType flag) { return static_cast<AtomType>(~static_cast<uint32_t>(flag)); }

// FIX: make consistent (naming, order, etc.)
namespace AtomTypeFlags {
inline bool has(AtomType flags, AtomType flag) {
  return (static_cast<uint32_t>(flags) & static_cast<uint32_t>(flag)) != 0;
}
inline bool has_any(AtomType flags, AtomType flag) { return (flags & flag) != AtomType::NONE; }
inline bool has_all(AtomType flags, AtomType flag) { return (flags & flag) == flag; }

inline AtomType remove(AtomType flags, AtomType flag) { return flags & ~flag; }

inline bool all(AtomType flags, AtomType toCheck) {
  return (static_cast<uint32_t>(flags) & static_cast<uint32_t>(toCheck)) == static_cast<uint32_t>(toCheck);
}

inline bool any(AtomType flags, AtomType toCheck) {
  return (static_cast<uint32_t>(flags) & static_cast<uint32_t>(toCheck)) != 0;
}

inline bool none(AtomType flags, AtomType toCheck) {
  return (static_cast<uint32_t>(flags) & static_cast<uint32_t>(toCheck)) == 0;
}

inline bool empty(AtomType flags) { return flags == AtomType::NONE; }

inline bool has_enum_as_string(AtomType flags, std::string flag_name) {
  return has(flags, string_to_atom_type(flag_name));
}

inline std::vector<AtomType> print_flags(AtomType flags) {
  std::vector<AtomType> result;
  for (int i = 0; i < 32; ++i) {
    AtomType flag = static_cast<AtomType>(1 << i);
    if (has(flags, flag)) {
      result.push_back(flag);
    }
  }
  return result;
}
} // namespace AtomTypeFlags

inline AtomType get_enum_as_string(std::string flag_name) { return string_to_atom_type(flag_name); }

inline std::string atom_type_to_string(AtomType type) {

  using namespace AtomTypeFlags;

  if (type == AtomType::NONE) return "None";

  std::string result;
  if (has(type, AtomType::HBOND_ACCEPTOR)) result += "HBOND_ACCEPTOR ";
  if (has(type, AtomType::HBOND_DONOR)) result += "HBOND_DONOR ";
  if (has(type, AtomType::WEAK_HBOND_ACCEPTOR)) result += "WEAK_HBOND_ACCEPTOR ";
  if (has(type, AtomType::WEAK_HBOND_DONOR)) result += "WEAK_HBOND_DONOR ";
  if (has(type, AtomType::POS_IONISABLE)) result += "POS_IONISABLE ";
  if (has(type, AtomType::NEG_IONISABLE)) result += "NEG_IONISABLE ";
  if (has(type, AtomType::CARBONYL_OXYGEN)) result += "CARBONYL_OXYGEN ";
  if (has(type, AtomType::CARBONYL_CARBON)) result += "CARBONYL_CARBON ";
  if (has(type, AtomType::AROMATIC)) result += "AROMATIC ";
  if (has(type, AtomType::HYDROPHOBIC)) result += "HYDROPHOBIC ";
  if (has(type, AtomType::XBOND_ACCEPTOR)) result += "XBOND_ACCEPTOR ";
  if (has(type, AtomType::XBOND_DONOR)) result += "XBOND_DONOR ";
  if (has(type, AtomType::IonicTypePartner)) result += "IONIC_TYPE_PARTNER ";
  if (has(type, AtomType::DativeBondPartner)) result += "DATIVE_PARTNER ";
  if (has(type, AtomType::TransitionMetal)) result += "TRANSITION_METAL ";
  if (has(type, AtomType::INVALID)) result += "UNKNOWN ";

  return result.empty() ? "Unknown" : result;
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
constexpr resTokenType operator"" _rt(char c) {
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
    if (entry == operator"" _rt(*values)) {
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

  resTokenType entry = lahuta::res_name_table(resname.c_str(), resname.length());

  // only standard amino acids
  if (static_cast<int>(entry) >= 20) {
    return AtomType::INVALID;
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
      return is_residue(entry, "N") ? 1033_at : 0_at;
    case encode_uint('N', 'E', '1'):
      return is_residue(entry, "W") ? 258_at : 0_at;
    case encode_uint('N', 'E', '2'):
      return is_residue(entry, "Q") ? 1031_at : is_residue(entry, "H") ? 1303_at : 0_at;
    case encode_uint('N', 'E'):
      return is_residue(entry, "R") ? 18_at : 0_at;
    case encode_uint('N', 'H', '1'):
    case encode_uint('N', 'H', '2'):
      return is_residue(entry, "R") ? 20_at : 0_at;
    case encode_uint('N', 'Z'):
      return is_residue(entry, "K") ? 18_at : 0_at;
    case encode_uint('O', 'D', '1'):
      return is_residue(entry, "N") ? 1031_at : is_residue(entry, "D") ? 1061_at : 0_at;
    case encode_uint('O', 'D', '2'):
      return is_residue(entry, "D") ? 1061_at : 0_at;
    case encode_uint('O', 'E', '1'):
      return is_residue(entry, "Q") ? 1031_at : is_residue(entry, "E") ? 1060_at : 0_at;
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

inline AtomType get_predef_aromatics(RDKit::Atom *at) {
  auto *info = static_cast<RDKit::AtomPDBResidueInfo *>(at->getMonomerInfo());
  resTokenType entry = res_name_table(info->getResidueName().c_str(), info->getResidueName().length());

  // Only consider standard residues
  if (static_cast<int>(entry) >= 20) {
    return AtomType::INVALID;
  }

  switch (entry) {
    case resTokenType::PHE:
    case resTokenType::TYR:
    case resTokenType::TRP:
    case resTokenType::HIS:
      break;
    default:
      return AtomType::NONE;
  }

  const std::string atom_name = info->getName();
  switch (encode_atom_name(atom_name)) {

    case encode_uint('C', 'D', '1'):
    case encode_uint('C', 'E', '2'):
      return !is_residue(entry, "H") ? AtomType::AROMATIC : AtomType::NONE;

    case encode_uint('C', 'D', '2'):
    case encode_uint('C', 'G'):
      return AtomType::AROMATIC;

    case encode_uint('C', 'E', '1'):
      return !is_residue(entry, "W") ? AtomType::AROMATIC : AtomType::NONE;

    case encode_uint('C', 'E', '3'):
    case encode_uint('C', 'Z', '2'):
    case encode_uint('C', 'Z', '3'):
    case encode_uint('C', 'H', '2'):
    case encode_uint('N', 'E', '1'):
      return is_residue(entry, "W") ? AtomType::AROMATIC : AtomType::NONE;

    case encode_uint('C', 'Z'):
      return is_residue(entry, "FY") ? AtomType::AROMATIC : AtomType::NONE;

    case encode_uint('N', 'D', '1'):
    case encode_uint('N', 'E', '2'):
      return is_residue(entry, "H") ? AtomType::AROMATIC : AtomType::NONE;

    default:
      return AtomType::NONE;
  }
}

constexpr std::pair<const char *, AtomType> AtomTypeSMARTS[] = {
    {"[$([nH]:@c(=O))]", 1029_at},
    {"[$([N;H2;v3;$(N-C(=O))])]", 1029_at},
    {"[$([n;H1;v3;!$([nH]cccc)])]", 1029_at},
    /*{"[#8,#9,$([#16;H0,H1;v2,v1]),$([N;v3;!$(N-*"*/
    /* "=!@[O,N,P,S]);!$(N-!@a);!$([NH]=!@*)]),$([nH0;+0])]",*/
    /* 1029_at},*/
    /**/
    {"[$(n:a:[nH])]", 2_at},
    {"[$([O;H0;$(O=C-[NH2])])]", 2_at},
    {"[$([O;H0;$(O=C([OH])-*)])]", 2_at},
    {"[N!H0v3,N!H0+v4,OH+0,SH+0,nH+0]", 2_at},

    {"[Cl,Br,I;X1;$([Cl,Br,I]-[#6])]", 10_at},

    {"[$([N;H2&+0][C;!$(C=*)]),$([N;H1&+0]([C;!$(C=*)])[C;!$(C=*)]),$([N;H0&+0]"
     "([C;!$(C=*)])([C;!$(C=*)])[C;!$(C=*)]);!$(N[a])]",
     16_at},
    {"NC(=N)", 16_at},
    {"[#7;+;!$([N+]-[O-])]", 16_at},
    {"[$([*+1,*+2,*+3]);!$([N+]-[O-])]", 16_at},

    {"[*-1,*-2]", 32_at},
    {"[$([OH,O-]-[C,S,N,P,Cl,Br,I]=O),$(O=[C,S,N,P,Cl,Br,I]-[OH,O-])]", 32_at},

    {"[$([OH0]=[CX3,c]);!$([OH0]=[CX3,c]-[OH,O-])]", 64_at},
    {"[$([CX3,c]=[OH0]);!$([CX3,c](=[OH0])-[OH,O-])]", 128_at},
    {"[#6+0!$(*~[#7,#8,F]),SH0+0v2,s+0,Cl+0,Br+0,I+0]", 512_at},

    // Require RingInfo initialization
    /*{"[n;R1]1[c;R1][n;R1][c;R1][c;R1]1", 16_at}, */
    /*{"[Xx]", 1_at},*/
};

inline std::vector<AtomType> match_atom_types(RDKit::ROMol &mol) {
  static std::array<RDKit::ROMol *, std::size(AtomTypeSMARTS)> patterns = [] {
    std::array<RDKit::ROMol *, std::size(AtomTypeSMARTS)> temp{};
    for (size_t i = 0; i < std::size(AtomTypeSMARTS); ++i) {
      temp[i] = RDKit::SmartsToMol(AtomTypeSMARTS[i].first);
    }
    return temp;
  }();

  RDKit::SubstructMatchParameters params;
  params.maxMatches = mol.getNumAtoms();

  std::vector<AtomType> types = {mol.getNumAtoms(), AtomType::NONE};
  for (size_t i = 0; i < std::size(AtomTypeSMARTS); ++i) {
    const auto &[smarts, atom_type] = AtomTypeSMARTS[i];
    RDKit::ROMol *pattern = patterns[i];

    SubStrMatches match_list;
    RDKit::SubstructMatch(mol, *pattern, match_list);

    for (const auto &match : match_list) {
      for (const auto &pair : match) {
        types[pair.second] |= atom_type;
      }
    }
  }

  return types;
}

} // namespace lahuta

#endif // ATOM_TYPES_HPP
