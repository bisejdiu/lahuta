#ifndef LAHUTA_TESTS_SASA_TEST_UTILS_HPP
#define LAHUTA_TESTS_SASA_TEST_UTILS_HPP

#include <array>
#include <ostream>

#include "analysis/sasa/sasa.hpp"

namespace lahuta::tests {

struct SasaMethodCase {
  const char *name = "";
  bool use_bitmask = false;
  bool use_simd    = false;
};

inline constexpr std::array<SasaMethodCase, 3> SasaMethodCases = {
    {
     {"standard", false, false},
     {"bitmask_no_simd", true, false},
     {"bitmask_simd", true, true},
     }
};

inline const std::array<SasaMethodCase, 3> &sasa_method_cases() { return SasaMethodCases; }

inline lahuta::analysis::SasaParams apply_sasa_method(const lahuta::analysis::SasaParams &base,
                                                      const SasaMethodCase &method) {
  auto params        = base;
  params.use_bitmask = method.use_bitmask;
  params.use_simd    = method.use_simd;
  return params;
}

inline std::ostream &operator<<(std::ostream &os, const SasaMethodCase &method) { return os << method.name; }

} // namespace lahuta::tests

#endif // LAHUTA_TESTS_SASA_TEST_UTILS_HPP
