#pragma once

#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <random>
#include <utility>
#include <vector>

// clang-format off
namespace lahuta::test_utils {
static_assert(std::numeric_limits<float>::is_iec559 && sizeof(float) == 4, "IEEE-754 32-bit float required");

// Deterministic mapping of a single 32-bit mt19937 draw to [0,1)
inline double u01(std::mt19937 &rng) noexcept {
  return static_cast<double>(rng()) * std::ldexp(1.0, -32);
}

inline double uniform(std::mt19937& rng, double lo, double hi) noexcept {
  const double t = u01(rng);
#if __cpp_lib_math_special_functions || (defined(__GNUC__) || defined(__clang__))
  return std::fma((hi - lo), t, lo);
#else
  return lo + (hi - lo) * t;
#endif
}

template <typename T>
inline std::vector<std::vector<T>> make_random_points(std::size_t n, std::mt19937 &rng, T lo, T hi) {
  std::vector<std::vector<T>> pts(n, std::vector<T>(3, T(0)));
  for (auto &p : pts) {
    p[0] = static_cast<T>(uniform(rng, static_cast<double>(lo), static_cast<double>(hi)));
    p[1] = static_cast<T>(uniform(rng, static_cast<double>(lo), static_cast<double>(hi)));
    p[2] = static_cast<T>(uniform(rng, static_cast<double>(lo), static_cast<double>(hi)));
  }
  return pts;
}

template <typename To, typename From>
inline std::vector<std::vector<To>> convert_vv(const std::vector<std::vector<From>> &src) {

  std::vector<std::vector<To>> dst;
  dst.reserve(src.size());
  for (const auto &row : src) {
    std::vector<To> r;
    r.reserve(3);
    r.push_back(static_cast<To>(row[0]));
    r.push_back(static_cast<To>(row[1]));
    r.push_back(static_cast<To>(row[2]));
    dst.push_back(std::move(r));
  }
  return dst;
}

inline long long ulp_distancef(float a, float b) {
  if (a == b) return 0;
  uint32_t ia, ib;
#if defined(__cpp_lib_bit_cast) && __cpp_lib_bit_cast >= 201806L
  ia = std::bit_cast<uint32_t>(a);
  ib = std::bit_cast<uint32_t>(b);
#else
  std::memcpy(&ia, &a, sizeof ia);
  std::memcpy(&ib, &b, sizeof ib);
#endif
  auto fix = [](uint32_t x) { return (x & 0x80000000u) ? (0x80000000u - x) : (x + 0x80000000u); };
  int64_t fa = static_cast<int64_t>(fix(ia));
  int64_t fb = static_cast<int64_t>(fix(ib));
  return (fa >= fb) ? (fa - fb) : (fb - fa);
}

[[nodiscard]] inline bool almost_equal_float(float a, float b, long long max_ulps, float abs_guard = 1e-6f) {
  if (a == b) return true;
  if (!std::isfinite(a) || !std::isfinite(b)) return a == b;
  if (std::fabs(static_cast<double>(a) - static_cast<double>(b)) <= abs_guard) return true;
  return ulp_distancef(a, b) <= max_ulps;
}

#define EXPECT_FLOAT_ULP_EQ(actual, expected, max_ulps, abs_guard, i, j)                                     \
  do {                                                                                                       \
    float _act = static_cast<float>(actual);                                                                 \
    float _exp = static_cast<float>(expected);                                                               \
    if (!lahuta::test_utils::almost_equal_float(_act, _exp, (max_ulps), (abs_guard))) {                      \
      auto _ulp = lahuta::test_utils::ulp_distancef(_act, _exp);                                             \
      ADD_FAILURE() << "Mismatch at (" << (i) << "," << (j) << ") actual=" << _act << " expect=" << _exp     \
                    << " ulp=" << _ulp;                                                                      \
    }                                                                                                        \
  } while (0)

} // namespace lahuta::test_utils
