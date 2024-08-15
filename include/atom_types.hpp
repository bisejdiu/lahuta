#include "bonds/token-gperf-generated.hpp"
#include <GraphMol/MonomerInfo.h>
#include <GraphMol/RWMol.h>
#include <cstdint>

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
};

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

inline bool hasFlag(AtomType flags, AtomType flag) {
  return (static_cast<uint32_t>(flags) & static_cast<uint32_t>(flag)) != 0;
}

inline bool hasAllFlags(AtomType flags, AtomType toCheck) {
  return (static_cast<uint32_t>(flags) & static_cast<uint32_t>(toCheck)) ==
         static_cast<uint32_t>(toCheck);
}

inline bool hasAnyFlag(AtomType flags, AtomType toCheck) {
  return (static_cast<uint32_t>(flags) & static_cast<uint32_t>(toCheck)) != 0;
}

inline bool hasNoFlag(AtomType flags, AtomType toCheck) {
  return (static_cast<uint32_t>(flags) & static_cast<uint32_t>(toCheck)) == 0;
}

inline bool getFlag(AtomType flags, AtomType flag) {
  return hasFlag(flags, flag);
}

inline void setFlag(AtomType &flags, AtomType flag) { flags |= flag; }

inline void clearFlag(AtomType &flags, AtomType flag) { flags &= ~flag; }

inline void toggleFlag(AtomType &flags, AtomType flag) { flags ^= flag; }

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
inline auto sfx = [](const std::string &name) {
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

  auto *info = static_cast<RDKit::AtomPDBResidueInfo *>(at->getMonomerInfo());

  auto entry = lahuta::res_name_table(info->getResidueName().c_str(),
                                      info->getResidueName().length());

  if (static_cast<int>(entry) >= 20) { // only standard amino acids handled
    return AtomType::NONE;
  }
  std::string name = info->getName();

  switch (sfx(name)) {
  case encode_uint('N'):
    if (info->getResidueName() == "PRO") {
      return AtomType::NONE;
    }
    return AtomType::HBOND_DONOR;
  case encode_uint('O'):
    return AtomType::HBOND_ACCEPTOR | AtomType::XBOND_ACCEPTOR |
           AtomType::WEAK_HBOND_ACCEPTOR | AtomType::CARBONYL_OXYGEN;
  case encode_uint('C'):
    return AtomType::CARBONYL_CARBON;
  case encode_uint('C', 'A'):
    return AtomType::WEAK_HBOND_DONOR;
  case encode_uint('C', 'B'):
    if (info->getResidueName() == "SER" || info->getResidueName() == "THR" ||
        info->getResidueName() == "TYR") {
      return AtomType::WEAK_HBOND_DONOR;
    }
    return AtomType::WEAK_HBOND_DONOR | AtomType::HYDROPHOBIC;
  case encode_uint('C', 'D'):
    if (info->getResidueName() == "LYS") {
      return AtomType::WEAK_HBOND_DONOR | AtomType::HYDROPHOBIC;
    }
    if (info->getResidueName() == "ARG" || info->getResidueName() == "PRO") {
      return AtomType::WEAK_HBOND_DONOR;
    }
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
  }
  case encode_uint('C', 'D', '2'):
    if (info->getResidueName() == "LYS") {
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
  case encode_uint('C', 'E', '2'):
    if (info->getResidueName() == "PHE" || info->getResidueName() == "TYR") {
      return AtomType::WEAK_HBOND_DONOR | AtomType::HYDROPHOBIC |
             AtomType::AROMATIC;
    }
    if (info->getResidueName() == "TRP") {
      return AtomType::AROMATIC;
    }
  case encode_uint('C', 'E', '3'):
    if (info->getResidueName() == "TRP") {
      return AtomType::HYDROPHOBIC | AtomType::AROMATIC;
    }
  case encode_uint('C', 'E'):
    if (info->getResidueName() == "MET") {
      return AtomType::WEAK_HBOND_DONOR | AtomType::HYDROPHOBIC;
    }
    if (info->getResidueName() == "LYS") {
      return AtomType::WEAK_HBOND_DONOR;
    }
  case encode_uint('C', 'G'):
    if (info->getResidueName() == "ARG" || info->getResidueName() == "GLN" ||
        info->getResidueName() == "GLU" || info->getResidueName() == "LEU" ||
        info->getResidueName() == "LYS" || info->getResidueName() == "MET" ||
        info->getResidueName() == "PRO") {
      return AtomType::WEAK_HBOND_DONOR | AtomType::HYDROPHOBIC;
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
  case encode_uint('C', 'G', '1'):
    if (info->getResidueName() == "ILE" || info->getResidueName() == "VAL") {
      return AtomType::WEAK_HBOND_DONOR | AtomType::HYDROPHOBIC;
    }
  case encode_uint('C', 'G', '2'):
    if (info->getResidueName() == "ILE" || info->getResidueName() == "VAL" ||
        info->getResidueName() == "THR") {
      return AtomType::WEAK_HBOND_DONOR | AtomType::HYDROPHOBIC;
    }
  case encode_uint('C', 'H', '2'):
    if (info->getResidueName() == "TRP") {
      return AtomType::WEAK_HBOND_DONOR | AtomType::HYDROPHOBIC |
             AtomType::AROMATIC;
    }
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
  case encode_uint('C', 'Z', '2'):
    if (info->getResidueName() == "TRP") {
      return AtomType::WEAK_HBOND_DONOR | AtomType::HYDROPHOBIC |
             AtomType::AROMATIC;
    }
  case encode_uint('C', 'Z', '3'):
    if (info->getResidueName() == "TRP") {
      return AtomType::WEAK_HBOND_DONOR | AtomType::HYDROPHOBIC |
             AtomType::AROMATIC;
    }
  case encode_uint('N', 'D', '1'):
    if (info->getResidueName() == "HIS") {
      return AtomType::HBOND_ACCEPTOR | AtomType::HBOND_DONOR |
             AtomType::XBOND_ACCEPTOR | AtomType::WEAK_HBOND_ACCEPTOR |
             AtomType::POS_IONISABLE | AtomType::AROMATIC;
    }
  case encode_uint('N', 'D', '2'):
    if (info->getResidueName() == "ASN") {
      return AtomType::HBOND_ACCEPTOR | AtomType::HBOND_DONOR |
             AtomType::XBOND_ACCEPTOR | AtomType::WEAK_HBOND_ACCEPTOR;
    }
  case encode_uint('N', 'E', '1'):
    if (info->getResidueName() == "TRP") {
      return AtomType::HBOND_DONOR | AtomType::AROMATIC;
    }
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
  case encode_uint('N', 'E'):
    if (info->getResidueName() == "ARG") {
      return AtomType::HBOND_DONOR | AtomType::POS_IONISABLE;
    }
  case encode_uint('N', 'H', '1'):
    if (info->getResidueName() == "ARG") {
      return AtomType::HBOND_DONOR | AtomType::POS_IONISABLE;
    }
  case encode_uint('N', 'H', '2'):
    if (info->getResidueName() == "ARG") {
      return AtomType::HBOND_DONOR | AtomType::POS_IONISABLE;
    }
  case encode_uint('N', 'Z'):
    if (info->getResidueName() == "LYS") {
      return AtomType::HBOND_DONOR | AtomType::POS_IONISABLE;
    }
  case encode_uint('O', 'D', '1'):
    if (info->getResidueName() == "ASN") {
      return AtomType::HBOND_ACCEPTOR | AtomType::HBOND_DONOR |
             AtomType::XBOND_ACCEPTOR | AtomType::WEAK_HBOND_ACCEPTOR;
    }
    if (info->getResidueName() == "ASP") {
      return AtomType::HBOND_ACCEPTOR | AtomType::XBOND_ACCEPTOR |
             AtomType::WEAK_HBOND_ACCEPTOR | AtomType::NEG_IONISABLE;
    }
  case encode_uint('O', 'D', '2'):
    if (info->getResidueName() == "ASP") {
      return AtomType::HBOND_ACCEPTOR | AtomType::XBOND_ACCEPTOR |
             AtomType::WEAK_HBOND_ACCEPTOR | AtomType::NEG_IONISABLE;
    }
  case encode_uint('O', 'E', '1'):
    if (info->getResidueName() == "GLN") {
      return AtomType::HBOND_ACCEPTOR | AtomType::HBOND_DONOR |
             AtomType::XBOND_ACCEPTOR | AtomType::WEAK_HBOND_ACCEPTOR;
    }
    if (info->getResidueName() == "GLU") {
      return AtomType::HBOND_ACCEPTOR | AtomType::XBOND_ACCEPTOR |
             AtomType::WEAK_HBOND_ACCEPTOR | AtomType::NEG_IONISABLE;
    }
  case encode_uint('O', 'E', '2'):
    if (info->getResidueName() == "GLU") {
      return AtomType::HBOND_ACCEPTOR | AtomType::XBOND_ACCEPTOR |
             AtomType::WEAK_HBOND_ACCEPTOR | AtomType::NEG_IONISABLE;
    }
  case encode_uint('O', 'G', '1'):
    if (info->getResidueName() == "THR") {
      return AtomType::HBOND_ACCEPTOR | AtomType::HBOND_DONOR |
             AtomType::XBOND_ACCEPTOR | AtomType::WEAK_HBOND_ACCEPTOR;
    }
  case encode_uint('O', 'G'):
    if (info->getResidueName() == "SER") {
      return AtomType::HBOND_ACCEPTOR | AtomType::HBOND_DONOR |
             AtomType::XBOND_ACCEPTOR | AtomType::WEAK_HBOND_ACCEPTOR;
    }
  case encode_uint('O', 'H'):
    if (info->getResidueName() == "TYR") {
      return AtomType::HBOND_ACCEPTOR | AtomType::HBOND_DONOR |
             AtomType::XBOND_ACCEPTOR | AtomType::WEAK_HBOND_ACCEPTOR;
    }
  case encode_uint('S', 'D'):
    if (info->getResidueName() == "MET") {
      return AtomType::HBOND_ACCEPTOR | AtomType::XBOND_ACCEPTOR |
             AtomType::WEAK_HBOND_ACCEPTOR | AtomType::HYDROPHOBIC;
    }
  case encode_uint('S', 'G'):
    if (info->getResidueName() == "CYS") {
      return AtomType::HBOND_ACCEPTOR | AtomType::HBOND_DONOR |
             AtomType::XBOND_ACCEPTOR | AtomType::WEAK_HBOND_ACCEPTOR;
    }
  default:
    return AtomType::NONE;
  }
}

inline std::string atom_type_to_string(AtomType type) {
  if (type == AtomType::NONE)
    return "None";

  std::string result;
  if (static_cast<uint32_t>(hasFlag(type, AtomType::HBOND_ACCEPTOR)))
    result += "HBOND_ACCEPTOR ";
  if (static_cast<uint32_t>(hasFlag(type, AtomType::HBOND_DONOR)))
    result += "HBOND_DONOR ";
  if (static_cast<uint32_t>(hasFlag(type, AtomType::WEAK_HBOND_ACCEPTOR)))
    result += "WEAK_HBOND_ACCEPTOR ";
  if (static_cast<uint32_t>(hasFlag(type, AtomType::WEAK_HBOND_DONOR)))
    result += "WEAK_HBOND_DONOR ";
  if (static_cast<uint32_t>(hasFlag(type, AtomType::POS_IONISABLE)))
    result += "POS_IONISABLE ";
  if (static_cast<uint32_t>(hasFlag(type, AtomType::NEG_IONISABLE)))
    result += "NEG_IONISABLE ";
  if (static_cast<uint32_t>(hasFlag(type, AtomType::CARBONYL_OXYGEN)))
    result += "CARBONYL_OXYGEN ";
  if (static_cast<uint32_t>(hasFlag(type, AtomType::CARBONYL_CARBON)))
    result += "CARBONYL_CARBON ";
  if (static_cast<uint32_t>(hasFlag(type, AtomType::AROMATIC)))
    result += "AROMATIC ";
  if (static_cast<uint32_t>(hasFlag(type, AtomType::HYDROPHOBIC)))
    result += "HYDROPHOBIC ";
  if (static_cast<uint32_t>(hasFlag(type, AtomType::XBOND_ACCEPTOR)))
    result += "XBOND_ACCEPTOR ";
  if (static_cast<uint32_t>(hasFlag(type, AtomType::XBOND_DONOR)))
    result += "XBOND_DONOR ";

  return result.empty() ? "Unknown" : result;
}
