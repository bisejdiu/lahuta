#ifndef LAHUTA_TOPOLOGY_FLAGS_HPP
#define LAHUTA_TOPOLOGY_FLAGS_HPP

#include <array>
#include <cstdint>
#include <type_traits>

// clang-format off
namespace lahuta {

enum class TopologyComputation : uint32_t {
  None              = 0,
  Neighbors         = 1 << 0,
  Bonds             = 1 << 1,
  NonStandardBonds  = 1 << 2,
  Residues          = 1 << 3,
  Rings             = 1 << 4,
  AtomTyping        = 1 << 5,

  // combinations
  Basic             = Neighbors | Bonds | NonStandardBonds,
  Standard          = Basic     | Residues,
  Extended          = Standard  | Rings,
  Complete          = Extended  | AtomTyping,
  All               = ~0u
};

/// generate bit flag array
template <typename E, unsigned N> constexpr std::array<E, N> make_bitmask_flags() {
  std::array<E, N> a{};
  for (unsigned i = 0; i < N; ++i) {
    a[i] = static_cast<E>(1u << i);
  }
  return a;
}

constexpr unsigned NUM_BASE_COMPUTATION_FLAGS = 6;

/// array of base flags
inline constexpr auto BASE_COMPUTATION_FLAGS = make_bitmask_flags<TopologyComputation, NUM_BASE_COMPUTATION_FLAGS>();

inline constexpr TopologyComputation operator|(TopologyComputation a, TopologyComputation b) {
  return static_cast<TopologyComputation>(
        static_cast<std::underlying_type_t<TopologyComputation>>(a)
      | static_cast<std::underlying_type_t<TopologyComputation>>(b));
}

inline constexpr TopologyComputation operator&(TopologyComputation a, TopologyComputation b) {
  return static_cast<TopologyComputation>(
        static_cast<std::underlying_type_t<TopologyComputation>>(a)
      & static_cast<std::underlying_type_t<TopologyComputation>>(b));
}

inline constexpr TopologyComputation operator~(TopologyComputation a) {
  return static_cast<TopologyComputation>(~static_cast<std::underlying_type_t<TopologyComputation>>(a));
}

inline constexpr bool has_flag(TopologyComputation flags, TopologyComputation flag) {
  return (flags & flag) == flag && flag != TopologyComputation::None;
}

/// Check if the value is a single base flag
inline constexpr bool is_base_flag(TopologyComputation comp) {
  for (auto flag : BASE_COMPUTATION_FLAGS) {
    if (comp == flag) return true;
  }
  return false;
}

} // namespace lahuta

#endif // LAHUTA_TOPOLOGY_FLAGS_HPP
