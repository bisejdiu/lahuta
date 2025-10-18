#ifndef LAHUTA_TYPING_FLAGS_HPP
#define LAHUTA_TYPING_FLAGS_HPP

#include <string>

#include "typing/types.hpp"

namespace lahuta {

using u32 = std::uint32_t;

// FIX: make consistent (naming, order, etc.)
namespace AtomTypeFlags {
inline bool has(AtomType flags, AtomType flag) {
  return (static_cast<u32>(flags) & static_cast<u32>(flag)) != 0;
}
inline bool has_any(AtomType flags, AtomType flag) { return (flags & flag) != AtomType::None; }
inline bool has_all(AtomType flags, AtomType flag) { return (flags & flag) == flag; }

inline AtomType remove(AtomType flags, AtomType flag) { 
  return static_cast<AtomType>(static_cast<uint32_t>(flags) & ~static_cast<uint32_t>(flag));
}

inline bool all(AtomType flags, AtomType flag) {
  return (static_cast<u32>(flags) & static_cast<u32>(flag)) == static_cast<u32>(flag);
}

inline bool any(AtomType flags, AtomType flag) {
  return (static_cast<u32>(flags) & static_cast<u32>(flag)) != 0;
}

inline bool none(AtomType flags, AtomType flag) {
  return (static_cast<u32>(flags) & static_cast<u32>(flag)) == 0;
}

inline bool empty(AtomType flags) { return flags == AtomType::None; }

inline std::vector<AtomType> split(AtomType flags) { // FIX: rename to get_components or similar
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


inline std::string atom_type_to_string(AtomType type) { // FIX: refactor
  using namespace AtomTypeFlags;

  if (type == AtomType::None) return "None";

  std::string result;
  if (has(type, AtomType::HbondAcceptor))     result += "HbondAcceptor ";
  if (has(type, AtomType::HbondDonor))        result += "HbondDonor ";
  if (has(type, AtomType::WeakHbondAcceptor)) result += "WeakHbondAcceptor ";
  if (has(type, AtomType::WeakHbondDonor))    result += "WeakHbondDonor ";
  if (has(type, AtomType::PositiveCharge))    result += "PositiveCharge ";
  if (has(type, AtomType::NegativeCharge))    result += "NegativeCharge ";
  if (has(type, AtomType::CarbonylOxygen))    result += "CarbonylOxygen ";
  if (has(type, AtomType::CarbonylCarbon))    result += "CarbonylCarbon ";
  if (has(type, AtomType::Aromatic))          result += "Aromatic ";
  if (has(type, AtomType::Hydrophobic))       result += "Hydrophobic ";
  if (has(type, AtomType::XBondAcceptor))     result += "XbondAcceptor ";
  if (has(type, AtomType::XbondDonor))        result += "XbondDonor ";
  if (has(type, AtomType::IonicTypePartner))  result += "IonicTypePartner ";
  if (has(type, AtomType::DativeBondPartner)) result += "DativeBondPartner ";
  if (has(type, AtomType::TransitionMetal))   result += "TransitionMetal ";
  if (has(type, AtomType::Invalid))           result += "Invalid ";

  if (!result.empty()) {
      result.pop_back();
      std::replace(result.begin(), result.end(), ' ', ',');
  }
  return result;
}

} // namespace lahuta

#endif // LAHUTA_TYPING_FLAGS_HPP
