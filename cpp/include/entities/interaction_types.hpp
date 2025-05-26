#ifndef LAHUTA_ENTITIES_INTERACTION_TYPES_HPP
#define LAHUTA_ENTITIES_INTERACTION_TYPES_HPP

#include <cstdint>

// clang-format off
namespace lahuta {

enum class Category : std::uint8_t {
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

enum class Flavor : std::uint8_t {
  Default,
  Parallel, // PiStacking “P”
  TShape    // PiStacking “T”
};

struct InteractionType {
  Category   category;
  Flavor     flavor;

  constexpr InteractionType(Category c = Category::None, Flavor f = Flavor::Default) noexcept
    : category(c), flavor(f) {}

  static const InteractionType None;
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

  constexpr operator std::uint8_t() const noexcept {
    return static_cast<std::uint8_t>(category) | (static_cast<std::uint8_t>(flavor) << 4);
  }

  constexpr bool operator==(const InteractionType& o) const noexcept {
    return category == o.category && flavor == o.flavor;
  }

  constexpr bool operator!=(const InteractionType& other) const noexcept {
    return !(*this == other);
  }
};

inline constexpr InteractionType InteractionType::None        {Category::None,        Flavor::Default};
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


} // namespace lahuta

#endif // LAHUTA_ENTITIES_INTERACTION_TYPES_HPP
