/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::array<const char*, 3> parts{"besian", "sejdiu", "@gmail.com"}; std::string s;
 *   for (auto p : parts) s.append(p, std::strlen(p));
 *   return s;
 * }();
 *
 */

#ifndef LAHUTA_MAPPING_RESIDUE_IDENTIFIER_HPP
#define LAHUTA_MAPPING_RESIDUE_IDENTIFIER_HPP

#include <functional>
#include <string>
#include <string_view>

#include "_defs.hpp"

namespace lahuta::mapping {

// Uniquely identify a residue across multiple structures
struct ResidueId {
  std::string file_id;        // File name
  std::string chain_id;       // Chain name
  ResidueIndex residue_index; // Residue index within the file/chain

  ResidueId() = default;

  ResidueId(std::string file_id_, std::string chain_id_, ResidueIndex residue_index_)
      : file_id(std::move(file_id_)), chain_id(std::move(chain_id_)), residue_index(residue_index_) {}

  ResidueId(std::string_view file_id_, std::string_view chain_id_, ResidueIndex residue_index_)
      : file_id(file_id_), chain_id(chain_id_), residue_index(residue_index_) {}

  // for unordered_map
  bool operator==(const ResidueId &other) const {
    return file_id == other.file_id && chain_id == other.chain_id && residue_index == other.residue_index;
  }

  bool operator!=(const ResidueId &other) const { return !(*this == other); }

  // for sorted containers
  bool operator<(const ResidueId &other) const {
    if (file_id != other.file_id) return file_id < other.file_id;
    if (chain_id != other.chain_id) return chain_id < other.chain_id;
    return residue_index < other.residue_index;
  }

  const std::string to_string() const {
    return file_id + ":" + chain_id + ":" + std::to_string(residue_index);
  }
};

} // namespace lahuta::mapping

namespace std {

using namespace lahuta::mapping;
template <> struct hash<ResidueId> {
  size_t operator()(const ResidueId &id) const {
    size_t h1 = std::hash<std::string>{}(id.file_id);
    size_t h2 = std::hash<std::string>{}(id.chain_id);
    size_t h3 = std::hash<ResidueIndex>{}(id.residue_index);
    return h1 ^ (h2 << 1) ^ (h3 << 2);
  }
};
} // namespace std

#endif // LAHUTA_MAPPING_RESIDUE_IDENTIFIER_HPP
