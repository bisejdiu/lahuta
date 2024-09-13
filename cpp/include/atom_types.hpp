#ifndef ATOM_TYPES_HPP
#define ATOM_TYPES_HPP

#include <cstdint>
#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/GraphMol/MonomerInfo.h>
#include <rdkit/GraphMol/RWMol.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/GraphMol/Substruct/SubstructMatch.h>

#include "bonds/token-gperf-generated.hpp"

namespace lahuta {

using SubStrMatches = std::vector<RDKit::MatchVectType>;

enum class AtomType : uint32_t {
  NONE = 0x0,
  HBOND_ACCEPTOR = 0x1,
  HBOND_DONOR = 0x2,
  WEAK_HBOND_ACCEPTOR = 0x4,
  WEAK_HBOND_DONOR = 0x8,
  POS_IONISABLE = 0x10,
  NEG_IONISABLE = 0x20,
  CARBONYL_OXYGEN = 0x40,
  CARBONYL_CARBON = 0x80,
  AROMATIC = 0x100,
  HYDROPHOBIC = 0x200,
  XBOND_ACCEPTOR = 0x400,
  XBOND_DONOR = 0x800,
  INVALID = 0x1000,
  // CYCLICAL = 0x1000,
  // INVALID = 0x2000,
  // HYDROGEN = 0x2000,
};

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
      // {"CYCLICAL", AtomType::CYCLICAL},
      {"INVALID", AtomType::INVALID}};

  auto it = stringToEnum.find(flag_name);
  if (it != stringToEnum.end()) {
    return it->second;
  }
  throw std::invalid_argument("Invalid AtomType flag name: " + flag_name);
}

inline AtomType operator|(AtomType lhs, AtomType rhs) {
  return static_cast<AtomType>(static_cast<uint32_t>(lhs) |
                               static_cast<uint32_t>(rhs));
}

inline AtomType &operator|=(AtomType &lhs, AtomType rhs) {
  lhs = lhs | rhs;
  return lhs;
}

inline AtomType operator&(AtomType lhs, AtomType rhs) {
  return static_cast<AtomType>(static_cast<uint32_t>(lhs) &
                               static_cast<uint32_t>(rhs));
}

inline AtomType &operator&=(AtomType &lhs, AtomType rhs) {
  lhs = lhs & rhs;
  return lhs;
}

inline AtomType operator^(AtomType lhs, AtomType rhs) {
  return static_cast<AtomType>(static_cast<uint32_t>(lhs) ^
                               static_cast<uint32_t>(rhs));
}

inline AtomType &operator^=(AtomType &lhs, AtomType rhs) {
  lhs = lhs ^ rhs;
  return lhs;
}

inline AtomType operator~(AtomType flag) {
  return static_cast<AtomType>(~static_cast<uint32_t>(flag));
}

namespace AtomTypeFlags {
inline bool has(AtomType flags, AtomType flag) {
  return (static_cast<uint32_t>(flags) & static_cast<uint32_t>(flag)) != 0;
}

inline bool has_enum_as_string(AtomType flags, std::string flag_name) {
  return has(flags, string_to_atom_type(flag_name));
}

inline AtomType remove(AtomType flags, AtomType flag) { return flags & ~flag; }

inline bool all(AtomType flags, AtomType toCheck) {
  return (static_cast<uint32_t>(flags) & static_cast<uint32_t>(toCheck)) ==
         static_cast<uint32_t>(toCheck);
}

inline bool any(AtomType flags, AtomType toCheck) {
  return (static_cast<uint32_t>(flags) & static_cast<uint32_t>(toCheck)) != 0;
}

inline bool none(AtomType flags, AtomType toCheck) {
  return (static_cast<uint32_t>(flags) & static_cast<uint32_t>(toCheck)) == 0;
}

inline bool empty(AtomType flags) { return flags == AtomType::NONE; }

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

inline AtomType get_enum_as_string(std::string flag_name) {
  return string_to_atom_type(flag_name);
}

inline std::string atom_type_to_string(AtomType type) {

  using namespace AtomTypeFlags;

  if (type == AtomType::NONE)
    return "None";

  std::string result;
  if (has(type, AtomType::HBOND_ACCEPTOR))
    result += "HBOND_ACCEPTOR ";
  if (has(type, AtomType::HBOND_DONOR))
    result += "HBOND_DONOR ";
  if (has(type, AtomType::WEAK_HBOND_ACCEPTOR))
    result += "WEAK_HBOND_ACCEPTOR ";
  if (has(type, AtomType::WEAK_HBOND_DONOR))
    result += "WEAK_HBOND_DONOR ";
  if (has(type, AtomType::POS_IONISABLE))
    result += "POS_IONISABLE ";
  if (has(type, AtomType::NEG_IONISABLE))
    result += "NEG_IONISABLE ";
  if (has(type, AtomType::CARBONYL_OXYGEN))
    result += "CARBONYL_OXYGEN ";
  if (has(type, AtomType::CARBONYL_CARBON))
    result += "CARBONYL_CARBON ";
  if (has(type, AtomType::AROMATIC))
    result += "AROMATIC ";
  if (has(type, AtomType::HYDROPHOBIC))
    result += "HYDROPHOBIC ";
  if (has(type, AtomType::XBOND_ACCEPTOR))
    result += "XBOND_ACCEPTOR ";
  if (has(type, AtomType::XBOND_DONOR))
    result += "XBOND_DONOR ";
  // if (has(type, AtomType::CYCLICAL))
  result += "CYCLICAL ";
  if (has(type, AtomType::INVALID))
    result += "UNKNOWN ";

  return result.empty() ? "Unknown" : result;
}

constexpr inline unsigned encode_uint(const char A, const char B = '\0',
                                      const char C = '\0') {
  unsigned result = static_cast<unsigned>(A) << 24;
  if (B != '\0') {
    result |= (static_cast<unsigned>(B) << 16);
    if (C != '\0') {
      result |= (static_cast<unsigned>(C) << 8);
    } else {
      result |= 0x1000000; // Flag for 2-character string
    }
  } else {
    result |= 0x2000000; // Flag for 1-character string
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

// FIX: This function does not need to check for non-standard residue; instead
// the caller should take over that responsibility.
static AtomType get_atom_type(RDKit::Atom *at) {
  // if (at->getAtomicNum() == 1) {
  //   return AtomType::HYDROGEN;
  // }

  auto *info = static_cast<RDKit::AtomPDBResidueInfo *>(at->getMonomerInfo());

  auto entry = lahuta::res_name_table(info->getResidueName().c_str(),
                                      info->getResidueName().length());

  if (static_cast<int>(entry) >= 20) { // only standard amino acids handled
    return AtomType::INVALID;
  }
  std::string name = info->getName();

  switch (encode_atom_name(name)) {
  case encode_uint('N'):
    if (info->getResidueName() == "PRO") {
      return AtomType::NONE;
      // return AtomType::CYCLICAL;
    }
    return AtomType::HBOND_DONOR;
  case encode_uint('O'):
    return AtomType::HBOND_ACCEPTOR | AtomType::XBOND_ACCEPTOR |
           AtomType::WEAK_HBOND_ACCEPTOR | AtomType::CARBONYL_OXYGEN;
  case encode_uint('C'):
    return AtomType::CARBONYL_CARBON;
  case encode_uint('C', 'A'):
    // if (info->getResidueName() == "PRO") {
    //   return AtomType::WEAK_HBOND_DONOR | AtomType::CYCLICAL;
    // }
    return AtomType::WEAK_HBOND_DONOR;
  case encode_uint('C', 'B'):
    if (info->getResidueName() == "SER" || info->getResidueName() == "THR" ||
        info->getResidueName() == "TYR") {
      return AtomType::WEAK_HBOND_DONOR;
    }
    // if (info->getResidueName() == "PRO") {
    //   return AtomType::WEAK_HBOND_DONOR | AtomType::HYDROPHOBIC |
    //   AtomType::CYCLICAL;
    // }
    return AtomType::WEAK_HBOND_DONOR | AtomType::HYDROPHOBIC;
  case encode_uint('C', 'D'):
    if (info->getResidueName() == "LYS") {
      return AtomType::WEAK_HBOND_DONOR | AtomType::HYDROPHOBIC;
    }
    if (info->getResidueName() == "ARG" || info->getResidueName() == "PRO") {
      return AtomType::WEAK_HBOND_DONOR;
    }
    // if (info->getResidueName() == "PRO") {
    //   return AtomType::CYCLICAL;
    // }
    return AtomType::NONE;
  case encode_uint('C', 'D', '1'): {
    if (info->getResidueName() == "TRP") {
      return AtomType::AROMATIC;
    }
    auto val = AtomType::WEAK_HBOND_DONOR | AtomType::HYDROPHOBIC;
    if (info->getResidueName() == "LEU" || info->getResidueName() == "ILE") {
      return val;
    }
    if (info->getResidueName() == "TYR" || info->getResidueName() == "PHE") {
      val |= AtomType::AROMATIC;
      return val;
    }
    return AtomType::NONE;
  }
  case encode_uint('C', 'D', '2'):
    if (info->getResidueName() == "LEU") {
      return AtomType::WEAK_HBOND_DONOR | AtomType::HYDROPHOBIC;
    }
    if (info->getResidueName() == "TYR" || info->getResidueName() == "PHE") {
      return AtomType::WEAK_HBOND_DONOR | AtomType::HYDROPHOBIC |
             AtomType::AROMATIC;
    }
    if (info->getResidueName() == "TRP") {
      return AtomType::HYDROPHOBIC | AtomType::AROMATIC;
    }
    if (info->getResidueName() == "HIS") {
      return AtomType::HBOND_ACCEPTOR | AtomType::HBOND_DONOR |
             AtomType::XBOND_ACCEPTOR | AtomType::WEAK_HBOND_ACCEPTOR |
             AtomType::POS_IONISABLE | AtomType::AROMATIC;
    }
    return AtomType::NONE;
  case encode_uint('C', 'E', '1'):
    if (info->getResidueName() == "HIS") {
      return AtomType::HBOND_ACCEPTOR | AtomType::HBOND_DONOR |
             AtomType::XBOND_ACCEPTOR | AtomType::WEAK_HBOND_ACCEPTOR |
             AtomType::POS_IONISABLE | AtomType::AROMATIC;
    }
    if (info->getResidueName() == "PHE" || info->getResidueName() == "TYR") {
      return AtomType::WEAK_HBOND_DONOR | AtomType::HYDROPHOBIC |
             AtomType::AROMATIC;
    }
    return AtomType::NONE;
  case encode_uint('C', 'E', '2'):
    if (info->getResidueName() == "PHE" || info->getResidueName() == "TYR") {
      return AtomType::WEAK_HBOND_DONOR | AtomType::HYDROPHOBIC |
             AtomType::AROMATIC;
    }
    if (info->getResidueName() == "TRP") {
      return AtomType::AROMATIC;
    }
    return AtomType::NONE;
  case encode_uint('C', 'E', '3'):
    if (info->getResidueName() == "TRP") {
      return AtomType::HYDROPHOBIC | AtomType::AROMATIC;
    }
    return AtomType::NONE;
  case encode_uint('C', 'E'):
    if (info->getResidueName() == "MET") {
      return AtomType::WEAK_HBOND_DONOR | AtomType::HYDROPHOBIC;
    }
    if (info->getResidueName() == "LYS") {
      return AtomType::WEAK_HBOND_DONOR;
    }
    return AtomType::NONE;
  case encode_uint('C', 'G'):
    if (info->getResidueName() == "ARG" || info->getResidueName() == "GLN" ||
        info->getResidueName() == "GLU" || info->getResidueName() == "LEU" ||
        info->getResidueName() == "LYS" || info->getResidueName() == "MET") {
      return AtomType::WEAK_HBOND_DONOR | AtomType::HYDROPHOBIC;
    }
    if (info->getResidueName() == "PRO") {
      return AtomType::WEAK_HBOND_DONOR |
             AtomType::HYDROPHOBIC; // | AtomType::CYCLICAL;
    }
    if (info->getResidueName() == "PHE") {
      return AtomType::WEAK_HBOND_DONOR | AtomType::HYDROPHOBIC |
             AtomType::AROMATIC;
    }
    if (info->getResidueName() == "TRP" || info->getResidueName() == "TYR") {
      return AtomType::HYDROPHOBIC | AtomType::AROMATIC;
    }
    if (info->getResidueName() == "HIS") {
      return AtomType::POS_IONISABLE | AtomType::AROMATIC;
    }
    return AtomType::NONE;
  case encode_uint('C', 'G', '1'):
    if (info->getResidueName() == "ILE" || info->getResidueName() == "VAL") {
      return AtomType::WEAK_HBOND_DONOR | AtomType::HYDROPHOBIC;
    }
    return AtomType::NONE;
  case encode_uint('C', 'G', '2'):
    if (info->getResidueName() == "ILE" || info->getResidueName() == "VAL" ||
        info->getResidueName() == "THR") {
      return AtomType::WEAK_HBOND_DONOR | AtomType::HYDROPHOBIC;
    }
    return AtomType::NONE;
  case encode_uint('C', 'H', '2'):
    if (info->getResidueName() == "TRP") {
      return AtomType::WEAK_HBOND_DONOR | AtomType::HYDROPHOBIC |
             AtomType::AROMATIC;
    }
    return AtomType::NONE;
  case encode_uint('C', 'Z'):
    if (info->getResidueName() == "ARG") {
      return AtomType::POS_IONISABLE;
    }
    if (info->getResidueName() == "PHE") {
      return AtomType::WEAK_HBOND_DONOR | AtomType::HYDROPHOBIC |
             AtomType::AROMATIC;
    }
    if (info->getResidueName() == "TYR") {
      return AtomType::AROMATIC;
    }
    return AtomType::NONE;
  case encode_uint('C', 'Z', '2'):
    if (info->getResidueName() == "TRP") {
      return AtomType::WEAK_HBOND_DONOR | AtomType::HYDROPHOBIC |
             AtomType::AROMATIC;
    }
    return AtomType::NONE;
  case encode_uint('C', 'Z', '3'):
    if (info->getResidueName() == "TRP") {
      return AtomType::WEAK_HBOND_DONOR | AtomType::HYDROPHOBIC |
             AtomType::AROMATIC;
    }
    return AtomType::NONE;
  case encode_uint('N', 'D', '1'):
    if (info->getResidueName() == "HIS") {
      return AtomType::HBOND_ACCEPTOR | AtomType::HBOND_DONOR |
             AtomType::XBOND_ACCEPTOR | AtomType::WEAK_HBOND_ACCEPTOR |
             AtomType::POS_IONISABLE | AtomType::AROMATIC;
    }
    return AtomType::NONE;
  case encode_uint('N', 'D', '2'):
    if (info->getResidueName() == "ASN") {
      return AtomType::HBOND_ACCEPTOR | AtomType::HBOND_DONOR |
             AtomType::XBOND_ACCEPTOR | AtomType::WEAK_HBOND_ACCEPTOR;
    }
    return AtomType::NONE;
  case encode_uint('N', 'E', '1'):
    if (info->getResidueName() == "TRP") {
      return AtomType::HBOND_DONOR | AtomType::AROMATIC;
    }
    return AtomType::NONE;
  case encode_uint('N', 'E', '2'):
    if (info->getResidueName() == "GLN") {
      return AtomType::HBOND_ACCEPTOR | AtomType::HBOND_DONOR |
             AtomType::XBOND_ACCEPTOR | AtomType::WEAK_HBOND_ACCEPTOR;
    }
    if (info->getResidueName() == "HIS") {
      return AtomType::HBOND_ACCEPTOR | AtomType::HBOND_DONOR |
             AtomType::XBOND_ACCEPTOR | AtomType::WEAK_HBOND_ACCEPTOR |
             AtomType::POS_IONISABLE | AtomType::AROMATIC;
    }
    return AtomType::NONE;
  case encode_uint('N', 'E'):
    if (info->getResidueName() == "ARG") {
      return AtomType::HBOND_DONOR | AtomType::POS_IONISABLE;
    }
    return AtomType::NONE;
  case encode_uint('N', 'H', '1'):
    if (info->getResidueName() == "ARG") {
      return AtomType::HBOND_DONOR | AtomType::POS_IONISABLE;
    }
    return AtomType::NONE;
  case encode_uint('N', 'H', '2'):
    if (info->getResidueName() == "ARG") {
      return AtomType::HBOND_DONOR | AtomType::POS_IONISABLE;
    }
    return AtomType::NONE;
  case encode_uint('N', 'Z'):
    if (info->getResidueName() == "LYS") {
      return AtomType::HBOND_DONOR | AtomType::POS_IONISABLE;
    }
    return AtomType::NONE;
  case encode_uint('O', 'D', '1'):
    if (info->getResidueName() == "ASN") {
      return AtomType::HBOND_ACCEPTOR | AtomType::HBOND_DONOR |
             AtomType::XBOND_ACCEPTOR | AtomType::WEAK_HBOND_ACCEPTOR;
    }
    if (info->getResidueName() == "ASP") {
      return AtomType::HBOND_ACCEPTOR | AtomType::XBOND_ACCEPTOR |
             AtomType::WEAK_HBOND_ACCEPTOR | AtomType::NEG_IONISABLE;
    }
    return AtomType::NONE;
  case encode_uint('O', 'D', '2'):
    if (info->getResidueName() == "ASP") {
      // if (at->getIdx() == 608) {
      //   std::cout << "O D 2: " << info->getResidueName() << " " <<
      //   info->getName() << std::endl;
      // }
      return AtomType::HBOND_ACCEPTOR | AtomType::XBOND_ACCEPTOR |
             AtomType::WEAK_HBOND_ACCEPTOR | AtomType::NEG_IONISABLE;
    }
    return AtomType::NONE;
  case encode_uint('O', 'E', '1'):
    if (info->getResidueName() == "GLN") {
      return AtomType::HBOND_ACCEPTOR | AtomType::HBOND_DONOR |
             AtomType::XBOND_ACCEPTOR | AtomType::WEAK_HBOND_ACCEPTOR;
    }
    if (info->getResidueName() == "GLU") {
      return AtomType::HBOND_ACCEPTOR | AtomType::XBOND_ACCEPTOR |
             AtomType::WEAK_HBOND_ACCEPTOR | AtomType::NEG_IONISABLE;
    }
    return AtomType::NONE;
  case encode_uint('O', 'E', '2'):
    if (info->getResidueName() == "GLU") {
      return AtomType::HBOND_ACCEPTOR | AtomType::XBOND_ACCEPTOR |
             AtomType::WEAK_HBOND_ACCEPTOR | AtomType::NEG_IONISABLE;
    }
    return AtomType::NONE;
  case encode_uint('O', 'G', '1'):
    if (info->getResidueName() == "THR") {
      return AtomType::HBOND_ACCEPTOR | AtomType::HBOND_DONOR |
             AtomType::XBOND_ACCEPTOR | AtomType::WEAK_HBOND_ACCEPTOR;
    }
    return AtomType::NONE;
  case encode_uint('O', 'G'):
    if (info->getResidueName() == "SER") {
      return AtomType::HBOND_ACCEPTOR | AtomType::HBOND_DONOR |
             AtomType::XBOND_ACCEPTOR | AtomType::WEAK_HBOND_ACCEPTOR;
    }
    return AtomType::NONE;
  case encode_uint('O', 'H'):
    if (info->getResidueName() == "TYR") {
      return AtomType::HBOND_ACCEPTOR | AtomType::HBOND_DONOR |
             AtomType::XBOND_ACCEPTOR | AtomType::WEAK_HBOND_ACCEPTOR;
    }
    return AtomType::NONE;
  case encode_uint('S', 'D'):
    if (info->getResidueName() == "MET") {
      return AtomType::HBOND_ACCEPTOR | AtomType::XBOND_ACCEPTOR |
             AtomType::WEAK_HBOND_ACCEPTOR | AtomType::HYDROPHOBIC;
    }
    return AtomType::NONE;
  case encode_uint('S', 'G'):
    if (info->getResidueName() == "CYS") {
      return AtomType::HBOND_ACCEPTOR | AtomType::HBOND_DONOR |
             AtomType::XBOND_ACCEPTOR | AtomType::WEAK_HBOND_ACCEPTOR;
    }
    return AtomType::NONE;
  case encode_uint('O', 'X', 'T'):
    return AtomType::HBOND_ACCEPTOR | AtomType::XBOND_ACCEPTOR |
           AtomType::WEAK_HBOND_ACCEPTOR;
  default:
    return AtomType::NONE;
  }
}

constexpr std::pair<const char *, AtomType> AtomTypeSMARTS[] = {
    {"[#8,#9,$([#16;H0,H1;v2,v1]),$([N;v3;!$(N-*=!@[O,N,P,S]);!$(N-!@a);!$([NH]"
     "=!@*)]),$([nH0;+0])]",
     AtomType::HBOND_ACCEPTOR},
    {"[$([nH]:@c(=O))]", AtomType::HBOND_ACCEPTOR},
    {"[$([n;H1;v3;!$([nH]cccc)])]", AtomType::HBOND_ACCEPTOR},
    {"[$([N;H2;v3;$(N-C(=O))])]", AtomType::HBOND_ACCEPTOR},

    {"[N!H0v3,N!H0+v4,OH+0,SH+0,nH+0]", AtomType::HBOND_DONOR},
    {"[$([O;H0;$(O=C([OH])-*)])]", AtomType::HBOND_DONOR},
    {"[$(n:a:[nH])]", AtomType::HBOND_DONOR},
    {"[$([O;H0;$(O=C-[NH2])])]", AtomType::HBOND_DONOR},

    // FIX: These seem to be very similar to weak hbond acceptor?
    // {"[#8,#9,$([#16;H0,H1;v2,v1]),$([N;v3;!$(N-*=!@[O,N,P,S]);!$(N-!@a);!$([NH]=!@*)]),$([nH0;+0])]",
    // AtomType::XBOND_ACCEPTOR},
    // {"[$([nH]:@c(=O))]", AtomType::XBOND_ACCEPTOR},
    // {"[$([n;H1;v3;!$([nH]cccc)])]", AtomType::XBOND_ACCEPTOR},
    // {"[$([N;H2;v3;$(N-C(=O))])]", AtomType::XBOND_ACCEPTOR},
    // {"[Xx]", AtomType::XBOND_ACCEPTOR},
    // {"[Cl,Br,I;X1;$([Cl,Br,I]-[#6])]", AtomType::XBOND_DONOR},

    {"[#8,#9,$([#16;H0,H1;v2,v1]),$([N;v3;!$(N-*=!@[O,N,P,S]);!$(N-!@a);!$([NH]"
     "=!@*)]),$([nH0;+0])]",
     AtomType::WEAK_HBOND_ACCEPTOR},
    {"[$([nH]:@c(=O))]", AtomType::WEAK_HBOND_ACCEPTOR},
    {"[$([n;H1;v3;!$([nH]cccc)])]", AtomType::WEAK_HBOND_ACCEPTOR},
    {"[$([N;H2;v3;$(N-C(=O))])]", AtomType::WEAK_HBOND_ACCEPTOR},
    {"[Cl,Br,I;X1;$([Cl,Br,I]-[#6])]", AtomType::WEAK_HBOND_ACCEPTOR},
    {"[#6!H0]", AtomType::WEAK_HBOND_DONOR},

    {"[$([N;H2&+0][C;!$(C=*)]),$([N;H1&+0]([C;!$(C=*)])[C;!$(C=*)]),$([N;H0&+0]"
     "([C;!$(C=*)])([C;!$(C=*)])[C;!$(C=*)]);!$(N[a])]",
     AtomType::POS_IONISABLE},
    {"[n;R1]1[c;R1][n;R1][c;R1][c;R1]1", AtomType::POS_IONISABLE},
    {"NC(=N)", AtomType::POS_IONISABLE},
    {"[#7;+;!$([N+]-[O-])]", AtomType::POS_IONISABLE},
    {"[$([*+1,*+2,*+3]);!$([N+]-[O-])]", AtomType::POS_IONISABLE},
    {"[Li,Be,Na,Mg,Al,K,Ca,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Ga,Rb,Sr,Y,Zr,Nb,Mo,Tc,"
     "Ru,Rh,Pd,Ag,Cd,In,Sn,Cs,Ba,La,Ce,Pr,Nd,Pm,Sm,Eu,Gd,Tb,Dy,Ho,Er,Tm,Yb,Lu,"
     "Hf,Ta,W,Re,Os,Ir,Pt,Au,Hg,Tl,Pb,Bi,Po,Fr,Ra,Ac,Th,Pa,U,Np,Pu,Am,Cm,Bk,"
     "Cf]",
     AtomType::POS_IONISABLE},

    {"[$([OH,O-]-[C,S,N,P,Cl,Br,I]=O),$(O=[C,S,N,P,Cl,Br,I]-[OH,O-])]",
     AtomType::NEG_IONISABLE},
    {"[*-1,*-2]", AtomType::NEG_IONISABLE},

    {"[#6+0!$(*~[#7,#8,F]),SH0+0v2,s+0,Cl+0,Br+0,I+0]", AtomType::HYDROPHOBIC},

    {"[$([OH0]=[CX3,c]);!$([OH0]=[CX3,c]-[OH,O-])]", AtomType::CARBONYL_OXYGEN},
    {"[$([CX3,c]=[OH0]);!$([CX3,c](=[OH0])-[OH,O-])]",
     AtomType::CARBONYL_CARBON},

    {"[a;r4,!R1&r3]1:[a;r4,!R1&r3]:[a;r4,!R1&r3]:[a;r4,!R1&r3]:1",
     AtomType::AROMATIC},
    {"[a;r5,!R1&r4,!R1&r3]1:[a;r5,!R1&r4,!R1&r3]:[a;r5,!R1&r4,!R1&r3]:[a;r5,!"
     "R1&r4,!R1&r3]:[a;r5,!R1&r4,!R1&r3]:1",
     AtomType::AROMATIC},
    {"[a;r6,!R1&r5,!R1&r4,!R1&r3]1:[a;r6,!R1&r5,!R1&r4,!R1&r3]:[a;r6,!R1&r5,!"
     "R1&r4,!R1&r3]:[a;r6,!R1&r5,!R1&r4,!R1&r3]:[a;r6,!R1&r5,!R1&r4,!R1&r3]:[a;"
     "r6,!R1&r5,!R1&r4,!R1&r3]:1",
     AtomType::AROMATIC},
    {"[a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]1:[a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:["
     "a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:[a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:[a;"
     "r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:[a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:[a;r7,"
     "!R1&r6,!R1&r5,!R1&r4,!R1&r3]:1",
     AtomType::AROMATIC},
    {"[a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]1:[a;r8,!R1&r7,!R1&r6,!R1&r5,!"
     "R1&r4,!R1&r3]:[a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:[a;r8,!R1&r7,!R1&"
     "r6,!R1&r5,!R1&r4,!R1&r3]:[a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:[a;r8,"
     "!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:[a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!"
     "R1&r3]:[a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:1",
     AtomType::AROMATIC},
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

    SubStrMatches matchList;
    RDKit::SubstructMatch(mol, *pattern, matchList);

    for (const auto &match : matchList) {
      for (const auto &pair : match) {
        types[pair.second] |= atom_type;
      }
      // auto *atom = mol.getAtomWithIdx(match[0].second);
      // if (atom->getAtomicNum() == 26) {
      //   auto *info =
      //   static_cast<RDKit::AtomPDBResidueInfo*>(atom->getMonomerInfo());
      //   std::cout << "-> Fe: " << info->getResidueName() << " "
      //             << info->getName() << " " << atom_type_to_string(atom_type)
      //             << std::endl;
      // }
      // types[match[0].second] |= atom_type;
    }
  }

  return types;
}

} // namespace lahuta

#endif // ATOM_TYPES_HPP
