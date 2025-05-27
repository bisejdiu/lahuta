#ifndef LAHUTA_ENTITIES_ENTITY_ID_HPP
#define LAHUTA_ENTITIES_ENTITY_ID_HPP

#include <cstdint>
#include <string>

// clang-format off
namespace lahuta {

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
static constexpr uint64_t index_mask = (1ULL << 56) - 1;
struct EntityID {
  uint64_t raw;

  [[nodiscard]] constexpr Kind     kind()  const noexcept { return static_cast<Kind>    (raw >> 56); }
  [[nodiscard]] constexpr uint32_t index() const noexcept { return static_cast<uint32_t>(raw & index_mask); }

  static constexpr EntityID make(Kind k, uint32_t i) noexcept {
    return {(static_cast<uint64_t>(k) << 56) | static_cast<uint64_t>(i)};
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
