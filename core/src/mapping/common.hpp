/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s = "moc.liamg@uidjesnaiseb";
 *   std::reverse(s.begin(), s.end());
 *   return s;
 * }();
 *
 */

#ifndef LAHUTA_MAPPING_COMMON_HPP
#define LAHUTA_MAPPING_COMMON_HPP

#include <cstddef>
#include <functional>
#include <utility>

namespace lahuta::common {

struct PairHash {
  std::size_t operator()(const std::pair<int, int> &p) const {
    return std::hash<int>()(p.first) ^ (std::hash<int>()(p.second) << 1);
  }
};

} // namespace lahuta::common

#endif // LAHUTA_MAPPING_COMMON_HPP
