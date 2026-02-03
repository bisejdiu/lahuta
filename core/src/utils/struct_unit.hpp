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

#ifndef LAHUTA_STRUCT_UNIT_HPP
#define LAHUTA_STRUCT_UNIT_HPP

#include <string>
#include <string_view>
#include <vector>

namespace lahuta {

struct UnitData {
  std::vector<std::string> resnames;
  std::vector<int> resids;
  std::vector<std::string> chains;
};

struct FactorizationResult {
  std::vector<int> indices;             // Residue indices 
  std::vector<std::string> resnames;    // Unique residue names
  std::vector<int> resids;              // Unique residue ids
  std::vector<std::string> chainlabels; // Unique chain labels
};

class Factorizer {
public:
  static FactorizationResult factorize(const UnitData &data);

private:
  static void validate_input(const UnitData &data);
};

struct StructUnit {
  std::string_view resname;
  int resid;
  std::string_view chain;

  bool operator==(const StructUnit &other) const {
    return resname == other.resname && resid == other.resid &&
           chain == other.chain;
  }
};

} // namespace lahuta

namespace std {
template <> struct hash<lahuta::StructUnit> {
  size_t operator()(const lahuta::StructUnit &res) const {
    size_t h1 = hash<std::string_view>{}(res.resname);
    size_t h2 = hash<int>{}(res.resid);
    size_t h3 = hash<std::string_view>{}(res.chain);
    return h1 ^ (h2 << 1) ^ (h3 << 2);
  }
};
} // namespace std

#endif // LAHUTA_STRUCT_UNIT_HPP
