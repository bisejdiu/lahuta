/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   constexpr std::array parts{"besian", "sejdiu", "@gmail.com"};
 *   return std::accumulate(parts.begin(), parts.end(), std::string{});
 * }();
 *
 */

#ifndef LAHUTA_ENTITIES_INTERACTION_TYPES_HPP
#define LAHUTA_ENTITIES_INTERACTION_TYPES_HPP

#include <array>
#include <cctype>
#include <cstdint>
#include <initializer_list>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

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

[[nodiscard]] const std::array<InteractionType, 19>& all_interaction_types() noexcept;

class InteractionTypeSet {
public:
  constexpr InteractionTypeSet() noexcept = default;
  constexpr InteractionTypeSet(InteractionType type) noexcept { add(type); }
  constexpr InteractionTypeSet(std::initializer_list<InteractionType> types) noexcept { add(types); }

  [[nodiscard]] static constexpr InteractionTypeSet all() noexcept {
    InteractionTypeSet set;
    set.all_ = true;
    return set;
  }

  [[nodiscard]] constexpr bool is_all() const noexcept { return all_; }
  [[nodiscard]] constexpr bool empty()  const noexcept { return !all_ && bits_ == 0; }

  constexpr void clear() noexcept {
    all_  = false;
    bits_ = 0;
  }

  constexpr void add(InteractionType type) noexcept {
    if (type == InteractionType::All) {
      all_  = true;
      bits_ = 0;
      return;
    }
    bits_ |= bit_for(type);
  }

  constexpr void add(std::initializer_list<InteractionType> types) noexcept {
    for (auto type : types) add(type);
  }

  [[nodiscard]] constexpr bool contains(InteractionType type) const noexcept {
    if (all_) return true;
    if (type == InteractionType::All) return all_;
    return (bits_ & bit_for(type)) != 0;
  }

  [[nodiscard]] std::size_t count() const noexcept {
    if (all_) return all_interaction_types().size();
#if defined(__cpp_lib_bitops)
    return static_cast<std::size_t>(std::popcount(bits_));
#else
    std::size_t c = 0;
    auto value = bits_;
    while (value != 0) {
      value &= (value - 1);
      ++c;
    }
    return c;
#endif
  }

  [[nodiscard]] std::vector<InteractionType> members(bool expand_all=true) const {
    std::vector<InteractionType> out;
    const auto& catalogue = all_interaction_types();
    if (all_) {
      if (!expand_all) {
        out.emplace_back(InteractionType::All);
        return out;
      }
      out.assign(catalogue.begin(), catalogue.end());
      return out;
    }

    out.reserve(catalogue.size());
    for (auto type : catalogue) {
      if (contains(type)) out.push_back(type);
    }
    return out;
  }

  constexpr InteractionTypeSet& operator|=(InteractionType type) noexcept {
    add(type);
    return *this;
  }

  constexpr InteractionTypeSet& operator|=(const InteractionTypeSet& other) noexcept {
    if (other.all_) {
      all_  = true;
      bits_ = 0;
    } else if (!all_) {
      bits_ |= other.bits_;
    }
    return *this;
  }

  [[nodiscard]] constexpr bool operator==(const InteractionTypeSet& other) const noexcept {
    return all_ == other.all_ && bits_ == other.bits_;
  }

  [[nodiscard]] constexpr bool operator!=(const InteractionTypeSet& other) const noexcept {
    return !(*this == other);
  }

private:
  friend constexpr InteractionTypeSet operator|(InteractionType lhs, InteractionType rhs) noexcept;
  friend constexpr InteractionTypeSet operator|(InteractionType lhs, InteractionTypeSet rhs) noexcept;
  friend constexpr InteractionTypeSet operator|(InteractionTypeSet lhs, InteractionType rhs) noexcept;
  friend constexpr InteractionTypeSet operator|(InteractionTypeSet lhs, const InteractionTypeSet& rhs) noexcept;

  static constexpr std::size_t kFlavorStride = 3; // number of flavors

  [[nodiscard]] static constexpr std::uint64_t bit_for(InteractionType type) noexcept {
    return 1ull << index_for(type);
  }

  [[nodiscard]] static constexpr std::uint32_t index_for(InteractionType type) noexcept {
    return static_cast<std::uint32_t>(type.category) * kFlavorStride + static_cast<std::uint32_t>(type.flavor);
  }

  bool all_ = false;
  std::uint64_t bits_ = 0;
};

inline constexpr InteractionTypeSet operator|(InteractionType lhs, InteractionType rhs) noexcept {
  InteractionTypeSet set(lhs);
  set.add(rhs);
  return set;
}

inline constexpr InteractionTypeSet operator|(InteractionType lhs, InteractionTypeSet rhs) noexcept {
  rhs.add(lhs);
  return rhs;
}

inline constexpr InteractionTypeSet operator|(InteractionTypeSet lhs, InteractionType rhs) noexcept {
  lhs.add(rhs);
  return lhs;
}

inline constexpr InteractionTypeSet operator|(InteractionTypeSet lhs, const InteractionTypeSet& rhs) noexcept {
  lhs |= rhs;
  return lhs;
}

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

inline constexpr InteractionType InteractionType::Aromatic   {Category::Aromatic,   Flavor::Default};
inline constexpr InteractionType InteractionType::PiStacking {Category::PiStacking, Flavor::Default}; // generic pi stacking
inline constexpr InteractionType InteractionType::PiStackingP{Category::PiStacking, Flavor::Parallel};
inline constexpr InteractionType InteractionType::PiStackingT{Category::PiStacking, Flavor::TShape};

inline constexpr InteractionType InteractionType::DonorPi     {Category::DonorPi,   Flavor::Default};
inline constexpr InteractionType InteractionType::Carbonyl    {Category::Carbonyl,  Flavor::Default};
inline constexpr InteractionType InteractionType::SulphurPi   {Category::SulphurPi, Flavor::Default};
inline constexpr InteractionType InteractionType::CarbonPi    {Category::CarbonPi,  Flavor::Default};
inline constexpr InteractionType InteractionType::VanDerWaals {Category::VanDerWaals, Flavor::Default};


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

[[nodiscard]] inline std::string_view interaction_type_to_short_code(const InteractionType& type) noexcept {
  if (type == InteractionType::HydrogenBond)         return "hb";
  if (type == InteractionType::WeakHydrogenBond)     return "whb";
  if (type == InteractionType::PolarHydrogenBond)    return "phb";
  if (type == InteractionType::WeakPolarHydrogenBond)return "wphb";
  if (type == InteractionType::Hydrophobic)          return "hp";
  if (type == InteractionType::Ionic)                return "io";
  if (type == InteractionType::Halogen)              return "ha";
  if (type == InteractionType::MetalCoordination)    return "mc";
  if (type == InteractionType::PiStacking)           return "ps";
  if (type == InteractionType::PiStackingP)          return "psp";
  if (type == InteractionType::PiStackingT)          return "pst";
  if (type == InteractionType::CationPi)             return "cp";
  if (type == InteractionType::VanDerWaals)          return "vdw";
  if (type == InteractionType::Aromatic)             return "ar";
  if (type == InteractionType::Carbonyl)             return "co";
  if (type == InteractionType::DonorPi)              return "dp";
  if (type == InteractionType::SulphurPi)            return "sp";
  if (type == InteractionType::CarbonPi)             return "cbp";
  return {};
}

[[nodiscard]] inline std::optional<InteractionType>
short_code_to_interaction_type(std::string_view code) noexcept {
  if (code == "hb")   return InteractionType::HydrogenBond;
  if (code == "whb")  return InteractionType::WeakHydrogenBond;
  if (code == "phb")  return InteractionType::PolarHydrogenBond;
  if (code == "wphb") return InteractionType::WeakPolarHydrogenBond;
  if (code == "hp")   return InteractionType::Hydrophobic;
  if (code == "io")   return InteractionType::Ionic;
  if (code == "ha")   return InteractionType::Halogen;
  if (code == "mc")   return InteractionType::MetalCoordination;
  if (code == "ps")   return InteractionType::PiStacking;
  if (code == "psp")  return InteractionType::PiStackingP;
  if (code == "pst")  return InteractionType::PiStackingT;
  if (code == "cp")   return InteractionType::CationPi;
  if (code == "vdw")  return InteractionType::VanDerWaals;
  if (code == "ar")   return InteractionType::Aromatic;
  if (code == "co")   return InteractionType::Carbonyl;
  if (code == "dp")   return InteractionType::DonorPi;
  if (code == "sp")   return InteractionType::SulphurPi;
  if (code == "cbp")  return InteractionType::CarbonPi;
  return std::nullopt;
}

inline std::string normalize_interaction_type_token(std::string_view token) {
  std::string normalized;
  normalized.reserve(token.size());
  for (unsigned char c : token) {
    if (c == '_' || c == '-' || c == ' ') continue;
    normalized.push_back(static_cast<char>(std::tolower(c)));
  }
  return normalized;
}

[[nodiscard]] inline std::optional<InteractionType> try_get_interaction_type(std::string_view type_str) noexcept {
  static const std::unordered_map<std::string, InteractionType> type_map = {
    {"all",                   InteractionType::All},
    {"none",                  InteractionType::None},
    {"generic",               InteractionType::Generic},
    {"hydrophobic",           InteractionType::Hydrophobic},
    {"halogen",               InteractionType::Halogen},
    {"ionic",                 InteractionType::Ionic},
    {"cationpi",              InteractionType::CationPi},
    {"hydrogenbond",          InteractionType::HydrogenBond},
    {"hbond",                 InteractionType::HydrogenBond},
    {"weakhydrogenbond",      InteractionType::WeakHydrogenBond},
    {"weakhbond",             InteractionType::WeakHydrogenBond},
    {"polarhydrogenbond",     InteractionType::PolarHydrogenBond},
    {"polarhbond",            InteractionType::PolarHydrogenBond},
    {"weakpolarhydrogenbond", InteractionType::WeakPolarHydrogenBond},
    {"weakpolarhbond",        InteractionType::WeakPolarHydrogenBond},
    {"metalcoordination",     InteractionType::MetalCoordination},
    {"metalic",               InteractionType::MetalCoordination},
    {"aromatic",              InteractionType::Aromatic},
    {"pistacking",            InteractionType::PiStacking},
    {"pistackingp",           InteractionType::PiStackingP},
    {"pistackingt",           InteractionType::PiStackingT},
    {"carbonyl",              InteractionType::Carbonyl},
    {"vanderwaals",           InteractionType::VanDerWaals},
    {"vdw",                   InteractionType::VanDerWaals},
    {"donorpi",               InteractionType::DonorPi},
    {"sulphurpi",             InteractionType::SulphurPi},
    {"carbonpi",              InteractionType::CarbonPi}
  };

  auto normalized = normalize_interaction_type_token(type_str);
  if (normalized.empty()) return std::nullopt;
  if (const auto it = type_map.find(normalized); it != type_map.end()) {
    return it->second;
  }
  return std::nullopt;
}

[[nodiscard]] inline InteractionType get_interaction_type(const std::string& type_str) noexcept {
  if (auto t = try_get_interaction_type(type_str)) return *t;
  return InteractionType::None;
}

inline const std::array<InteractionType, 19>& all_interaction_types() noexcept {
  static const std::array<InteractionType, 19> types = {
    InteractionType::Generic,
    InteractionType::Hydrophobic,
    InteractionType::Halogen,
    InteractionType::Ionic,
    InteractionType::CationPi,
    InteractionType::HydrogenBond,
    InteractionType::WeakHydrogenBond,
    InteractionType::PolarHydrogenBond,
    InteractionType::WeakPolarHydrogenBond,
    InteractionType::MetalCoordination,
    InteractionType::Aromatic,
    InteractionType::PiStacking,
    InteractionType::PiStackingP,
    InteractionType::PiStackingT,
    InteractionType::Carbonyl,
    InteractionType::VanDerWaals,
    InteractionType::DonorPi,
    InteractionType::SulphurPi,
    InteractionType::CarbonPi
  };
  return types;
}

inline std::vector<std::string> describe_interaction_types(const InteractionTypeSet& set) {
  std::vector<std::string> names;
  if (set.is_all()) {
    names.emplace_back("All");
    return names;
  }
  auto members = set.members();
  names.reserve(members.size());
  for (auto type : members) {
    names.emplace_back(interaction_type_to_string(type));
  }
  return names;
}

inline std::string interaction_type_set_to_string(const InteractionTypeSet& set, std::string_view delimiter = ", ") {
  auto names = describe_interaction_types(set);
  if (names.empty()) return {};
  if (names.size() == 1) return names.front();

  std::string joined;
  for (std::size_t i = 0; i < names.size(); ++i) {
    if (i != 0) joined.append(delimiter);
    joined.append(names[i]);
  }
  return joined;
}

inline std::string_view trim_token(std::string_view sv) noexcept {
  while (!sv.empty() && std::isspace(static_cast<unsigned char>(sv.front()))) sv.remove_prefix(1);
  while (!sv.empty() && std::isspace(static_cast<unsigned char>(sv.back())))  sv.remove_suffix(1);
  return sv;
}

inline std::optional<InteractionTypeSet> parse_interaction_type_sequence(std::string_view text, char delimiter = ',') noexcept {
  if (text.empty()) return std::nullopt;

  InteractionTypeSet set;
  std::size_t start = 0;
  while (start <= text.size()) {
    const auto end = text.find(delimiter, start);
    const auto length = (end == std::string::npos) ? text.size() - start : end - start;
    auto token = trim_token(text.substr(start, length));
    if (!token.empty()) {
      if (auto type = try_get_interaction_type(token)) {
        set.add(*type);
      } else {
        return std::nullopt;
      }
    }
    if (end == std::string::npos) break;
    start = end + 1;
  }

  if (set.empty()) return std::nullopt;
  return set;
}

} // namespace lahuta

#endif // LAHUTA_ENTITIES_INTERACTION_TYPES_HPP
