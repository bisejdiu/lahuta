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
