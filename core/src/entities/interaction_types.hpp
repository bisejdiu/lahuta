#ifndef LAHUTA_ENTITIES_INTERACTION_TYPES_HPP
#define LAHUTA_ENTITIES_INTERACTION_TYPES_HPP

#include <cstdint>
#include <string>
#include <unordered_map>

// clang-format off
namespace lahuta {

enum class Category : std::uint16_t {
  None,
  Generic,
  Hydrophobic,
  Halogen,
  HydrogenBond,
  WeakHydrogenBond,
  PolarHydrogenBond,
  WeakPolarHydrogenBond,
  Aromatic,
  Ionic,
  MetalCoordination,
  CationPi,
  PiStacking,
  Carbonyl,
  VanDerWaals,
  DonorPi,
  SulphurPi,
  CarbonPi
};

enum class Flavor : std::uint16_t {
  Default,
  Parallel, // PiStacking P
  TShape    // PiStacking T
};

struct InteractionType {
  Category   category;
  Flavor     flavor;

  constexpr InteractionType(Category c = Category::None, Flavor f = Flavor::Default) noexcept
    : category(c), flavor(f) {}

  static const InteractionType None; // internal no-hit marker
  static const InteractionType All;
  static const InteractionType Generic;
  static const InteractionType Hydrophobic;
  static const InteractionType Halogen;
  static const InteractionType Ionic;
  static const InteractionType CationPi;
  static const InteractionType HydrogenBond;
  static const InteractionType WeakHydrogenBond;
  static const InteractionType PolarHydrogenBond;
  static const InteractionType WeakPolarHydrogenBond;
  static const InteractionType MetalCoordination;
  static const InteractionType Aromatic;
  static const InteractionType PiStacking;
  static const InteractionType PiStackingP;
  static const InteractionType PiStackingT;

  static const InteractionType Carbonyl;
  static const InteractionType VanDerWaals;
  static const InteractionType DonorPi;
  static const InteractionType SulphurPi;
  static const InteractionType CarbonPi;

  // 32-bit packed code: [ flavor:16 | category:16 ]
  constexpr operator std::uint32_t() const noexcept {
    return static_cast<std::uint32_t>(category) | (static_cast<std::uint32_t>(flavor) << 16);
  }

  constexpr bool operator==(const InteractionType& o) const noexcept {
    return category == o.category && flavor == o.flavor;
  }

  constexpr bool operator!=(const InteractionType& other) const noexcept {
    return !(*this == other);
  }
};

inline constexpr InteractionType InteractionType::None        {Category::None,        Flavor::Default};
inline constexpr InteractionType InteractionType::All         {Category::None,        Flavor::Parallel};
inline constexpr InteractionType InteractionType::Generic     {Category::Generic,     Flavor::Default};
inline constexpr InteractionType InteractionType::Hydrophobic {Category::Hydrophobic, Flavor::Default};
inline constexpr InteractionType InteractionType::Halogen     {Category::Halogen,     Flavor::Default};
inline constexpr InteractionType InteractionType::Ionic       {Category::Ionic,       Flavor::Default};
inline constexpr InteractionType InteractionType::CationPi    {Category::CationPi,    Flavor::Default};

inline constexpr InteractionType InteractionType::HydrogenBond     {Category::HydrogenBond,      Flavor::Default};
inline constexpr InteractionType InteractionType::WeakHydrogenBond {Category::WeakHydrogenBond,  Flavor::Default};
inline constexpr InteractionType InteractionType::PolarHydrogenBond     {Category::PolarHydrogenBond,      Flavor::Default};
inline constexpr InteractionType InteractionType::WeakPolarHydrogenBond {Category::WeakPolarHydrogenBond,  Flavor::Default};
inline constexpr InteractionType InteractionType::MetalCoordination{Category::MetalCoordination, Flavor::Default};

inline constexpr InteractionType InteractionType::Aromatic     {Category::Aromatic, Flavor::Default};
inline constexpr InteractionType InteractionType::PiStacking {Category::PiStacking, Flavor::Default}; // generic pi stacking
inline constexpr InteractionType InteractionType::PiStackingP{Category::PiStacking, Flavor::Parallel};
inline constexpr InteractionType InteractionType::PiStackingT{Category::PiStacking, Flavor::TShape};

inline constexpr InteractionType InteractionType::Carbonyl     {Category::Carbonyl, Flavor::Default};
inline constexpr InteractionType InteractionType::VanDerWaals {Category::VanDerWaals, Flavor::Default};
inline constexpr InteractionType InteractionType::DonorPi     {Category::DonorPi, Flavor::Default};
inline constexpr InteractionType InteractionType::SulphurPi   {Category::SulphurPi, Flavor::Default};
inline constexpr InteractionType InteractionType::CarbonPi    {Category::CarbonPi, Flavor::Default};


[[nodiscard]] inline std::string interaction_type_to_string(const InteractionType& type) noexcept {
  if (type == InteractionType::All) return "All";
  switch (type.category) {
    case Category::None:              return "None";
    case Category::Generic:           return "Generic";
    case Category::Hydrophobic:       return "Hydrophobic";
    case Category::Halogen:           return "Halogen";
    case Category::HydrogenBond:          return "HydrogenBond";
    case Category::WeakHydrogenBond:      return "WeakHydrogenBond";
    case Category::PolarHydrogenBond:     return "PolarHydrogenBond";
    case Category::WeakPolarHydrogenBond: return "WeakPolarHydrogenBond";
    case Category::MetalCoordination:     return "MetalCoordination";
    case Category::Aromatic:      return "Aromatic";
    case Category::Ionic:         return "Ionic";
    case Category::CationPi:      return "CationPi";
    case Category::PiStacking:
      if (type.flavor == Flavor::Parallel) return "PiStackingP";
      if (type.flavor == Flavor::TShape)   return "PiStackingT";
      return "PiStacking";
    case Category::Carbonyl:    return "Carbonyl";
    case Category::VanDerWaals: return "VanDerWaals";
    case Category::DonorPi:     return "DonorPi";
    case Category::SulphurPi:   return "SulphurPi";
    case Category::CarbonPi:    return "CarbonPi";
    default: return "Unknown";
  }
}

[[nodiscard]] inline InteractionType get_interaction_type(const std::string& type_str) noexcept {
  static const std::unordered_map<std::string, InteractionType> type_map = {
    {"all",              InteractionType::All},
    {"hbond",            InteractionType::HydrogenBond},
    {"weak_hbond",       InteractionType::WeakHydrogenBond},
    {"polar_hbond",      InteractionType::PolarHydrogenBond},
    {"weak_polar_hbond", InteractionType::WeakPolarHydrogenBond},
    {"hydrophobic",      InteractionType::Hydrophobic},
    {"ionic",            InteractionType::Ionic},
    {"halogen",          InteractionType::Halogen},
    {"metalic",          InteractionType::MetalCoordination},
    {"cationpi",         InteractionType::CationPi},
    {"pistacking",       InteractionType::PiStacking},
    {"aromatic",         InteractionType::Aromatic},
    {"carbonyl",         InteractionType::Carbonyl},
    {"vdw",              InteractionType::VanDerWaals},
    {"donor_pi",         InteractionType::DonorPi},
    {"sulphur_pi",       InteractionType::SulphurPi},
    {"carbon_pi",        InteractionType::CarbonPi}
  };

  if (const auto it = type_map.find(type_str); it != type_map.end()) {
      return it->second;
  }
  return InteractionType::None;
}

} // namespace lahuta

#endif // LAHUTA_ENTITIES_INTERACTION_TYPES_HPP
