#ifndef LAHUTA_HASH_FNV1A_HPP
#define LAHUTA_HASH_FNV1A_HPP

#include <cstddef>
#include <cstdint>
#include <string>

// clang-format off
namespace lahuta {

constexpr uint64_t fnv1a_64(const char *s, std::size_t len) noexcept {
  uint64_t hash = 14695981039346656037ULL;
  for (std::size_t i = 0; i < len; ++i) {
    hash ^= static_cast<uint64_t>(static_cast<unsigned char>(s[i]));
    hash *= 1099511628211ULL;
  }
  return hash;
}

inline uint64_t fnv1a_64(const std::string& s) noexcept {
  return fnv1a_64(s.data(), s.size());
}

} // namespace lahuta

#endif // LAHUTA_HASH_FNV1A_HPP
