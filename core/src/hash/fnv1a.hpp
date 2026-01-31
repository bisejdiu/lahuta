#ifndef LAHUTA_HASH_FNV1A_HPP
#define LAHUTA_HASH_FNV1A_HPP

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>

namespace lahuta {

// FNV-1a 32-bit
constexpr uint32_t fnv1a_32(const char *s, std::size_t len) noexcept {
  uint32_t hash = 2166136261u;
  for (std::size_t i = 0; i < len; ++i) {
    hash ^= static_cast<uint32_t>(static_cast<unsigned char>(s[i]));
    hash *= 16777619u;
  }
  return hash;
}

inline uint32_t fnv1a_32(std::string_view s) noexcept { return fnv1a_32(s.data(), s.size()); }

// FNV-1a 64-bit
constexpr uint64_t fnv1a_64(const char *s, std::size_t len) noexcept {
  uint64_t hash = 14695981039346656037ULL;
  for (std::size_t i = 0; i < len; ++i) {
    hash ^= static_cast<uint64_t>(static_cast<unsigned char>(s[i]));
    hash *= 1099511628211ULL;
  }
  return hash;
}

inline uint64_t fnv1a_64(const std::string &s) noexcept { return fnv1a_64(s.data(), s.size()); }

} // namespace lahuta

#endif // LAHUTA_HASH_FNV1A_HPP
