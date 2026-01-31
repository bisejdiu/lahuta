#ifndef LAHUTA_ANALYSIS_SASA_BITMASK_LUT_HPP
#define LAHUTA_ANALYSIS_SASA_BITMASK_LUT_HPP

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "analysis/sasa/sphere.hpp"
#include "analysis/sasa/types.hpp"

namespace lahuta::analysis {

inline constexpr int BitmaskDirResolution = 24;
inline constexpr int BitmaskAngleBins     = 256;
static_assert(BitmaskAngleBins > 1, "BitmaskAngleBins must be > 1");

struct BitmaskLut {
  SpherePointCount n_points = 0;
  int dir_resolution        = 0;
  int angle_bins            = 0;
  std::size_t words         = 0;
  std::vector<double> dir_x;
  std::vector<double> dir_y;
  std::vector<double> dir_z;
  std::vector<std::uint64_t> masks;
  std::uint64_t last_word_mask = 0;
};

namespace detail {

[[nodiscard]] inline int clamp_int(int value, int lo, int hi) noexcept {
  if (value < lo) return lo;
  if (value > hi) return hi;
  return value;
}

struct UnitVec {
  double x = 0.0;
  double y = 0.0;
  double z = 1.0;
};

[[nodiscard]] inline UnitVec octa_decode(double u, double v) noexcept {
  double x = u;
  double y = v;
  double z = 1.0 - std::abs(u) - std::abs(v);
  if (z < 0.0) {
    const double x_sign = (x >= 0.0) ? 1.0 : -1.0;
    const double y_sign = (y >= 0.0) ? 1.0 : -1.0;
    const double x_abs  = std::abs(x);
    const double y_abs  = std::abs(y);

    x = (1.0 - y_abs) * x_sign;
    y = (1.0 - x_abs) * y_sign;
  }
  const double norm = std::sqrt(x * x + y * y + z * z);
  if (norm <= 0.0) return {};

  return {x / norm, y / norm, z / norm};
}

inline void octa_encode(double x, double y, double z, double &u, double &v) noexcept {
  const double inv_l1 = 1.0 / (std::abs(x) + std::abs(y) + std::abs(z));

  double px = x * inv_l1;
  double py = y * inv_l1;
  double pz = z * inv_l1;

  if (pz < 0.0) {
    const double px_sign = (px >= 0.0) ? 1.0 : -1.0;
    const double py_sign = (py >= 0.0) ? 1.0 : -1.0;
    const double px_abs  = std::abs(px);
    const double py_abs  = std::abs(py);

    px = (1.0 - py_abs) * px_sign;
    py = (1.0 - px_abs) * py_sign;
  }
  u = px;
  v = py;
}

} // namespace detail

// Compute direction bin index from a unit vector using octahedral encoding.
[[nodiscard]] inline int dir_bin_from_unit(double x, double y, double z, int resolution) noexcept {
  double u = 0.0;
  double v = 0.0;
  detail::octa_encode(x, y, z, u, v);
  const double fx = (u + 1.0) * 0.5 * static_cast<double>(resolution);
  const double fy = (v + 1.0) * 0.5 * static_cast<double>(resolution);
  const int ix    = detail::clamp_int(static_cast<int>(fx), 0, resolution - 1);
  const int iy    = detail::clamp_int(static_cast<int>(fy), 0, resolution - 1);
  return iy * resolution + ix;
}

[[nodiscard]] inline int angle_bin_from_cos(double cos_value, int bins) noexcept {
  const double scaled = (cos_value + 1.0) * 0.5 * static_cast<double>(bins - 1);
  int idx             = static_cast<int>(std::lround(scaled));
  return detail::clamp_int(idx, 0, bins - 1);
}

// Population count for 64-bit integers
[[nodiscard]] inline std::size_t popcount64(std::uint64_t v) noexcept {
#if defined(__GNUG__) || defined(__clang__)
  return static_cast<std::size_t>(__builtin_popcountll(v));
#else
  std::size_t count = 0;
  while (v) {
    v &= (v - 1);
    ++count;
  }
  return count;
#endif
}

[[nodiscard]] inline bool supports_bitmask(SpherePointCount n_points) noexcept {
  return n_points == 64 || n_points == 128 || n_points == 256;
}

[[nodiscard]] inline BitmaskLut build_bitmask_lut(SpherePointCount n_points) {
  BitmaskLut lut;
  lut.n_points       = n_points;
  lut.dir_resolution = BitmaskDirResolution;
  lut.angle_bins     = BitmaskAngleBins;
  lut.words          = (n_points + 63) / 64;

  const std::size_t remainder = n_points % 64;
  lut.last_word_mask          = (remainder == 0) ? ~std::uint64_t{0} : ((std::uint64_t{1} << remainder) - 1);

  const int dir_count = lut.dir_resolution * lut.dir_resolution;
  lut.dir_x.resize(dir_count);
  lut.dir_y.resize(dir_count);
  lut.dir_z.resize(dir_count);

  for (int iy = 0; iy < lut.dir_resolution; ++iy) {
    for (int ix = 0; ix < lut.dir_resolution; ++ix) {
      const double u = (static_cast<double>(ix) + 0.5) / static_cast<double>(lut.dir_resolution) * 2.0 - 1.0;
      const double v = (static_cast<double>(iy) + 0.5) / static_cast<double>(lut.dir_resolution) * 2.0 - 1.0;
      const auto dir = detail::octa_decode(u, v);
      const int idx  = iy * lut.dir_resolution + ix;
      lut.dir_x[idx] = dir.x;
      lut.dir_y[idx] = dir.y;
      lut.dir_z[idx] = dir.z;
    }
  }

  const auto sphere            = generate_sphere(n_points);
  const std::size_t mask_count = static_cast<std::size_t>(dir_count) * lut.angle_bins * lut.words;
  lut.masks.assign(mask_count, 0);

  for (int dir_idx = 0; dir_idx < dir_count; ++dir_idx) {
    const double dx = lut.dir_x[dir_idx];
    const double dy = lut.dir_y[dir_idx];
    const double dz = lut.dir_z[dir_idx];
    for (int angle_idx = 0; angle_idx < lut.angle_bins; ++angle_idx) {
      const double cos_value = -1.0 + (2.0 * static_cast<double>(angle_idx)) /
                                          static_cast<double>(lut.angle_bins - 1);
      std::uint64_t *mask = lut.masks.data() +
                            (static_cast<std::size_t>(dir_idx) * lut.angle_bins + angle_idx) * lut.words;
      for (std::size_t k = 0; k < n_points; ++k) {
        const double dot = sphere.x[k] * dx + sphere.y[k] * dy + sphere.z[k] * dz;
        if (dot <= cos_value) {
          const std::size_t word  = k / 64;
          const std::size_t bit   = k % 64;
          mask[word]             |= (std::uint64_t{1} << bit);
        }
      }
      if (lut.words > 0) {
        mask[lut.words - 1] &= lut.last_word_mask;
      }
    }
  }

  return lut;
}

[[nodiscard]] inline const BitmaskLut *get_bitmask_lut(SpherePointCount n_points) noexcept {
  if (n_points == 64) {
    static const BitmaskLut lut = build_bitmask_lut(64);
    return &lut;
  }
  if (n_points == 128) {
    static const BitmaskLut lut = build_bitmask_lut(128);
    return &lut;
  }
  if (n_points == 256) {
    static const BitmaskLut lut = build_bitmask_lut(256);
    return &lut;
  }
  return nullptr;
}

} // namespace lahuta::analysis

#endif // LAHUTA_ANALYSIS_SASA_BITMASK_LUT_HPP
