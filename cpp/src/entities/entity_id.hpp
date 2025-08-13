#ifndef LAHUTA_ENTITIES_ENTITY_ID_HPP
#define LAHUTA_ENTITIES_ENTITY_ID_HPP

#include <cstdint>
#include <string>

// clang-format off
namespace lahuta {

using u64 = std::uint64_t;
using u32 = std::uint32_t;

enum class Kind : uint8_t { Atom = 0, Ring = 1, Group = 2 }; // order encodes priority

inline const char *kind_to_string(Kind kind) {
  switch (kind) {
    case Kind::Atom:  return "Atom";
    case Kind::Ring:  return "Ring";
    case Kind::Group: return "Group";
    default: throw std::invalid_argument("Invalid Kind value");
  }
}

// The upper 8 bits store the Kind, and the lower 56 bits store the index.
static constexpr u64 index_mask = (1ULL << 56) - 1;
struct EntityID {
  u64 raw;

  [[nodiscard]] constexpr Kind  kind() const noexcept { return static_cast<Kind>(raw >> 56); }
  [[nodiscard]] constexpr u32  index() const noexcept { return static_cast<u32> (raw & index_mask); }

  static constexpr EntityID make(Kind k, u32 i) noexcept {
    return {(static_cast<u64>(k) << 56) | static_cast<u64>(i)};
  }

  std::string to_string() const {
    std::string result = kind_to_string(kind());
    result += "#";
    result += std::to_string(index());
    return result;
  }

  friend constexpr bool operator==(EntityID a, EntityID b) noexcept { return a.raw == b.raw; }
  friend constexpr bool operator!=(EntityID a, EntityID b) noexcept { return a.raw != b.raw; }
  friend constexpr bool operator< (EntityID a, EntityID b) noexcept { return a.raw <  b.raw; }
};

} // namespace lahuta

#endif // LAHUTA_ENTITIES_ENTITY_ID_HPP
