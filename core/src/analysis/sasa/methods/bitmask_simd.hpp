/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s(22, '\0'); char* ptr = s.data();
 *   for (char c : std::string_view{"besian"}) *ptr++ = c;
 *   for (char c : std::string_view{"sejdiu"}) *ptr++ = c;
 *   *ptr++ = '@';
 *   for (char c : std::string_view{"gmail.com"}) *ptr++ = c;
 *   return s;
 * }();
 *
 */

#ifndef LAHUTA_ANALYSIS_SASA_METHODS_BITMASK_SIMD_HPP
#define LAHUTA_ANALYSIS_SASA_METHODS_BITMASK_SIMD_HPP

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <stdexcept>

#include <simde/x86/avx2.h>

#include "analysis/sasa/bitmask_lut.hpp"
#include "analysis/sasa/method.hpp"
#include "utils/math_constants.hpp"

namespace lahuta::analysis::simd {

using lahuta::Pi;

//
// I tried fast_rsqrt (using the Quake III "magic number" approximation).
// However, benchmarking showed:
//   - Fast rsqrt SIMD:     16-30% speedup, 97-99% exact match,  0.22%   deviation
//   - High-precision SIMD: 16-30% speedup, 99.98% exact match, <0.0001% deviation
//
// The performance difference was negligible because the main SIMD gains come
// from 4x neighbor batching and reduced loop overhead. I am choosing high-precision
// so users can use SIMD without reservation.   - Besian, Jan 2026
//

[[nodiscard]] inline double fast_rsqrt_scalar(double x) noexcept { return 1.0 / std::sqrt(x); }

[[nodiscard]] inline simde__m256d fast_rsqrt_avx2(simde__m256d x) noexcept {
  const simde__m256d one = simde_mm256_set1_pd(1.0);
  return simde_mm256_div_pd(one, simde_mm256_sqrt_pd(x));
}

// Vectorized octahedral encoding for 4 unit vectors.
// Given normalized direction vectors (nx, ny, nz), computes octahedral (u, v) coordinates.
inline void octa_encode_avx2(simde__m256d nx, simde__m256d ny, simde__m256d nz, simde__m256d &u,
                             simde__m256d &v) noexcept {
  const simde__m256d sign_mask = simde_mm256_set1_pd(-0.0);
  const simde__m256d one       = simde_mm256_set1_pd(1.0);
  const simde__m256d zero      = simde_mm256_setzero_pd();

  simde__m256d abs_nx = simde_mm256_andnot_pd(sign_mask, nx);
  simde__m256d abs_ny = simde_mm256_andnot_pd(sign_mask, ny);
  simde__m256d abs_nz = simde_mm256_andnot_pd(sign_mask, nz);

  simde__m256d l1_norm = simde_mm256_add_pd(abs_nx, simde_mm256_add_pd(abs_ny, abs_nz));
  simde__m256d inv_l1  = simde_mm256_div_pd(one, l1_norm);

  simde__m256d px = simde_mm256_mul_pd(nx, inv_l1);
  simde__m256d py = simde_mm256_mul_pd(ny, inv_l1);
  simde__m256d pz = simde_mm256_mul_pd(nz, inv_l1);

  simde__m256d abs_px = simde_mm256_andnot_pd(sign_mask, px);
  simde__m256d abs_py = simde_mm256_andnot_pd(sign_mask, py);

  simde__m256d one_minus_abs_py = simde_mm256_sub_pd(one, abs_py);
  simde__m256d one_minus_abs_px = simde_mm256_sub_pd(one, abs_px);

  simde__m256d px_sign = simde_mm256_and_pd(px, sign_mask);
  simde__m256d py_sign = simde_mm256_and_pd(py, sign_mask);

  simde__m256d flip_x = simde_mm256_or_pd(one_minus_abs_py, px_sign);
  simde__m256d flip_y = simde_mm256_or_pd(one_minus_abs_px, py_sign);

  simde__m256d neg_mask = simde_mm256_cmp_pd(pz, zero, SIMDE_CMP_LT_OQ);

  u = simde_mm256_blendv_pd(px, flip_x, neg_mask);
  v = simde_mm256_blendv_pd(py, flip_y, neg_mask);
}

// Compute direction bin indices for 4 vectors. Returns indices in [0, resolution^2).
inline void dir_bin_avx2(simde__m256d u, simde__m256d v, int resolution, int *indices) noexcept {
  const double half_res           = static_cast<double>(resolution) * 0.5;
  const simde__m256d half_res_vec = simde_mm256_set1_pd(half_res);

  simde__m256d fx = simde_mm256_add_pd(simde_mm256_mul_pd(u, half_res_vec), half_res_vec);
  simde__m256d fy = simde_mm256_add_pd(simde_mm256_mul_pd(v, half_res_vec), half_res_vec);

  simde__m128i ix_lo = simde_mm256_cvttpd_epi32(fx);
  simde__m128i iy_lo = simde_mm256_cvttpd_epi32(fy);

  alignas(16) int ix_arr[4];
  alignas(16) int iy_arr[4];
  simde_mm_storeu_si128(reinterpret_cast<simde__m128i *>(ix_arr), ix_lo);
  simde_mm_storeu_si128(reinterpret_cast<simde__m128i *>(iy_arr), iy_lo);

  for (int i = 0; i < 4; ++i) {
    int ix     = ix_arr[i];
    int iy     = iy_arr[i];
    ix         = (ix < 0) ? 0 : ((ix >= resolution) ? resolution - 1 : ix);
    iy         = (iy < 0) ? 0 : ((iy >= resolution) ? resolution - 1 : iy);
    indices[i] = iy * resolution + ix;
  }
}

// Compute angle bin indices for 4 cos_threshold values.
inline void angle_bin_avx2(simde__m256d cos_threshold, int bins, int *indices) noexcept {
  const double scale           = 0.5 * static_cast<double>(bins - 1);
  const simde__m256d one       = simde_mm256_set1_pd(1.0);
  const simde__m256d scale_vec = simde_mm256_set1_pd(scale);
  const simde__m256d half      = simde_mm256_set1_pd(0.5);

  simde__m256d cos_plus_one = simde_mm256_add_pd(cos_threshold, one);
  simde__m256d scaled       = simde_mm256_add_pd(simde_mm256_mul_pd(cos_plus_one, scale_vec), half);

  simde__m128i idx_lo = simde_mm256_cvttpd_epi32(scaled);

  alignas(16) int idx_arr[4];
  simde_mm_storeu_si128(reinterpret_cast<simde__m128i *>(idx_arr), idx_lo);

  for (int i = 0; i < 4; ++i) {
    int idx    = idx_arr[i];
    idx        = (idx < 0) ? 0 : ((idx >= bins) ? bins - 1 : idx);
    indices[i] = idx;
  }
}

// Batched kernel for 64-point bitmask (1 word).
template <typename CoordAccessor>
[[nodiscard]] double compute_bitmask_64_batched(const BitmaskLut &lut,
                                                const ComputeContext<CoordAccessor> &ctx) {
  const double r_i        = ctx.radius();
  const double r_i_sq     = ctx.radius_sq();
  const std::size_t start = ctx.neighbor_start();
  const std::size_t end   = ctx.neighbor_end();

  if (start == end) return 4.0 * Pi * r_i_sq;

  const double inv_two_r = 1.0 / (2.0 * r_i);
  const double cx        = ctx.coords.x(ctx.atom_index);
  const double cy        = ctx.coords.y(ctx.atom_index);
  const double cz        = ctx.coords.z(ctx.atom_index);

  std::uint64_t visibility = lut.last_word_mask;

  const simde__m256d cx_vec        = simde_mm256_set1_pd(cx);
  const simde__m256d cy_vec        = simde_mm256_set1_pd(cy);
  const simde__m256d cz_vec        = simde_mm256_set1_pd(cz);
  const simde__m256d r_i_sq_vec    = simde_mm256_set1_pd(r_i_sq);
  const simde__m256d inv_two_r_vec = simde_mm256_set1_pd(inv_two_r);
  const simde__m256d neg_one       = simde_mm256_set1_pd(-1.0);
  const simde__m256d pos_one       = simde_mm256_set1_pd(1.0);

  std::size_t idx = start;

  while (idx + 4 <= end) {
    alignas(32) double nbr_cx[4];
    alignas(32) double nbr_cy[4];
    alignas(32) double nbr_cz[4];
    alignas(32) double nbr_thresh_sq[4];

    for (int b = 0; b < 4; ++b) {
      const auto &nbr     = ctx.neighbors.items[idx + b];
      const std::size_t j = nbr.index;
      nbr_cx[b]           = ctx.coords.x(j);
      nbr_cy[b]           = ctx.coords.y(j);
      nbr_cz[b]           = ctx.coords.z(j);
      nbr_thresh_sq[b]    = nbr.threshold_sq;
    }

    simde__m256d nbr_cx_vec    = simde_mm256_load_pd(nbr_cx);
    simde__m256d nbr_cy_vec    = simde_mm256_load_pd(nbr_cy);
    simde__m256d nbr_cz_vec    = simde_mm256_load_pd(nbr_cz);
    simde__m256d thresh_sq_vec = simde_mm256_load_pd(nbr_thresh_sq);

    simde__m256d vx    = simde_mm256_sub_pd(cx_vec, nbr_cx_vec);
    simde__m256d vy    = simde_mm256_sub_pd(cy_vec, nbr_cy_vec);
    simde__m256d vz    = simde_mm256_sub_pd(cz_vec, nbr_cz_vec);
    simde__m256d vx_sq = simde_mm256_mul_pd(vx, vx);
    simde__m256d vy_sq = simde_mm256_mul_pd(vy, vy);
    simde__m256d vz_sq = simde_mm256_mul_pd(vz, vz);

    simde__m256d v_mag_sq  = simde_mm256_add_pd(vx_sq, simde_mm256_add_pd(vy_sq, vz_sq));
    simde__m256d numerator = simde_mm256_sub_pd(simde_mm256_sub_pd(thresh_sq_vec, v_mag_sq), r_i_sq_vec);
    simde__m256d limit     = simde_mm256_mul_pd(numerator, inv_two_r_vec);

    simde__m256d inv_dist      = fast_rsqrt_avx2(v_mag_sq);
    simde__m256d cos_threshold = simde_mm256_mul_pd(limit, inv_dist);

    simde__m256d full_occlude_mask = simde_mm256_cmp_pd(cos_threshold, pos_one, SIMDE_CMP_GE_OQ);
    if (simde_mm256_movemask_pd(full_occlude_mask) != 0) {
      alignas(32) double cos_arr[4];
      simde_mm256_store_pd(cos_arr, cos_threshold);
      for (int b = 0; b < 4; ++b) {
        if (cos_arr[b] >= 1.0) {
          return 0.0;
        }
      }
    }

    simde__m256d valid_mask = simde_mm256_cmp_pd(cos_threshold, neg_one, SIMDE_CMP_GT_OQ);
    int valid_bits          = simde_mm256_movemask_pd(valid_mask);

    if (valid_bits != 0) {
      simde__m256d nx = simde_mm256_mul_pd(vx, inv_dist);
      simde__m256d ny = simde_mm256_mul_pd(vy, inv_dist);
      simde__m256d nz = simde_mm256_mul_pd(vz, inv_dist);

      simde__m256d u, v;
      octa_encode_avx2(nx, ny, nz, u, v);

      alignas(16) int dir_indices[4];
      alignas(16) int angle_indices[4];
      dir_bin_avx2(u, v, lut.dir_resolution, dir_indices);
      angle_bin_avx2(cos_threshold, lut.angle_bins, angle_indices);

      std::uint64_t combined_mask = 0;
      for (int b = 0; b < 4; ++b) {
        if (valid_bits & (1 << b)) {
          const std::size_t mask_offset = static_cast<std::size_t>(dir_indices[b]) * lut.angle_bins +
                                          angle_indices[b];
          combined_mask |= lut.masks[mask_offset];
        }
      }

      visibility &= ~combined_mask;
      if (visibility == 0) return 0.0;
    }

    idx += 4;
  }

  while (idx < end) {
    const auto &nbr       = ctx.neighbors.items[idx];
    const std::size_t j   = nbr.index;
    const double vx       = cx - ctx.coords.x(j);
    const double vy       = cy - ctx.coords.y(j);
    const double vz       = cz - ctx.coords.z(j);
    const double v_mag_sq = vx * vx + vy * vy + vz * vz;

    if (v_mag_sq == 0.0) {
      throw std::invalid_argument("SASA invalid geometry: coincident atom centers.");
    }

    const double inv_dist      = fast_rsqrt_scalar(v_mag_sq);
    const double limit         = (nbr.threshold_sq - v_mag_sq - r_i_sq) * inv_two_r;
    const double cos_threshold = limit * inv_dist;

    if (cos_threshold >= 1.0) return 0.0;

    if (cos_threshold > -1.0) {
      const double nx                = vx * inv_dist;
      const double ny                = vy * inv_dist;
      const double nz                = vz * inv_dist;
      const int dir_idx              = dir_bin_from_unit(nx, ny, nz, lut.dir_resolution);
      const int angle_idx            = angle_bin_from_cos(cos_threshold, lut.angle_bins);
      const std::size_t mask_offset  = static_cast<std::size_t>(dir_idx) * lut.angle_bins + angle_idx;
      visibility                    &= ~lut.masks[mask_offset];

      if (visibility == 0) return 0.0;
    }
    ++idx;
  }

  const std::size_t accessible = popcount64(visibility);
  return (4.0 * Pi * r_i_sq / 64.0) * static_cast<double>(accessible);
}

// Batched kernel for 128-point bitmask (2 words).
template <typename CoordAccessor>
[[nodiscard]] double compute_bitmask_128_batched(const BitmaskLut &lut,
                                                 const ComputeContext<CoordAccessor> &ctx) {
  const double r_i        = ctx.radius();
  const double r_i_sq     = ctx.radius_sq();
  const std::size_t start = ctx.neighbor_start();
  const std::size_t end   = ctx.neighbor_end();

  if (start == end) return 4.0 * Pi * r_i_sq;

  const double inv_two_r = 1.0 / (2.0 * r_i);
  const double cx        = ctx.coords.x(ctx.atom_index);
  const double cy        = ctx.coords.y(ctx.atom_index);
  const double cz        = ctx.coords.z(ctx.atom_index);

  alignas(16) std::uint64_t vis_arr[2] = {~std::uint64_t{0}, lut.last_word_mask};
  simde__m128i visibility              = simde_mm_load_si128(reinterpret_cast<const simde__m128i *>(vis_arr));

  const simde__m256d cx_vec        = simde_mm256_set1_pd(cx);
  const simde__m256d cy_vec        = simde_mm256_set1_pd(cy);
  const simde__m256d cz_vec        = simde_mm256_set1_pd(cz);
  const simde__m256d r_i_sq_vec    = simde_mm256_set1_pd(r_i_sq);
  const simde__m256d inv_two_r_vec = simde_mm256_set1_pd(inv_two_r);
  const simde__m256d neg_one       = simde_mm256_set1_pd(-1.0);
  const simde__m256d pos_one       = simde_mm256_set1_pd(1.0);

  std::size_t idx = start;

  while (idx + 4 <= end) {
    alignas(32) double nbr_cx[4];
    alignas(32) double nbr_cy[4];
    alignas(32) double nbr_cz[4];
    alignas(32) double nbr_thresh_sq[4];

    for (int b = 0; b < 4; ++b) {
      const auto &nbr     = ctx.neighbors.items[idx + b];
      const std::size_t j = nbr.index;
      nbr_cx[b]           = ctx.coords.x(j);
      nbr_cy[b]           = ctx.coords.y(j);
      nbr_cz[b]           = ctx.coords.z(j);
      nbr_thresh_sq[b]    = nbr.threshold_sq;
    }

    simde__m256d nbr_cx_vec    = simde_mm256_load_pd(nbr_cx);
    simde__m256d nbr_cy_vec    = simde_mm256_load_pd(nbr_cy);
    simde__m256d nbr_cz_vec    = simde_mm256_load_pd(nbr_cz);
    simde__m256d thresh_sq_vec = simde_mm256_load_pd(nbr_thresh_sq);

    simde__m256d vx    = simde_mm256_sub_pd(cx_vec, nbr_cx_vec);
    simde__m256d vy    = simde_mm256_sub_pd(cy_vec, nbr_cy_vec);
    simde__m256d vz    = simde_mm256_sub_pd(cz_vec, nbr_cz_vec);
    simde__m256d vx_sq = simde_mm256_mul_pd(vx, vx);
    simde__m256d vy_sq = simde_mm256_mul_pd(vy, vy);
    simde__m256d vz_sq = simde_mm256_mul_pd(vz, vz);

    simde__m256d v_mag_sq  = simde_mm256_add_pd(vx_sq, simde_mm256_add_pd(vy_sq, vz_sq));
    simde__m256d numerator = simde_mm256_sub_pd(simde_mm256_sub_pd(thresh_sq_vec, v_mag_sq), r_i_sq_vec);
    simde__m256d limit     = simde_mm256_mul_pd(numerator, inv_two_r_vec);

    simde__m256d inv_dist      = fast_rsqrt_avx2(v_mag_sq);
    simde__m256d cos_threshold = simde_mm256_mul_pd(limit, inv_dist);

    simde__m256d full_occlude_mask = simde_mm256_cmp_pd(cos_threshold, pos_one, SIMDE_CMP_GE_OQ);
    if (simde_mm256_movemask_pd(full_occlude_mask) != 0) {
      alignas(32) double cos_arr[4];
      simde_mm256_store_pd(cos_arr, cos_threshold);
      for (int b = 0; b < 4; ++b) {
        if (cos_arr[b] >= 1.0) {
          return 0.0;
        }
      }
    }

    simde__m256d valid_mask = simde_mm256_cmp_pd(cos_threshold, neg_one, SIMDE_CMP_GT_OQ);
    int valid_bits          = simde_mm256_movemask_pd(valid_mask);

    if (valid_bits != 0) {
      simde__m256d nx = simde_mm256_mul_pd(vx, inv_dist);
      simde__m256d ny = simde_mm256_mul_pd(vy, inv_dist);
      simde__m256d nz = simde_mm256_mul_pd(vz, inv_dist);

      simde__m256d u, v;
      octa_encode_avx2(nx, ny, nz, u, v);

      alignas(16) int dir_indices[4];
      alignas(16) int angle_indices[4];
      dir_bin_avx2(u, v, lut.dir_resolution, dir_indices);
      angle_bin_avx2(cos_threshold, lut.angle_bins, angle_indices);

      simde__m128i combined_mask = simde_mm_setzero_si128();
      for (int b = 0; b < 4; ++b) {
        if (valid_bits & (1 << b)) {
          const std::size_t mask_offset = (static_cast<std::size_t>(dir_indices[b]) * lut.angle_bins +
                                           angle_indices[b]) *
                                          2;
          const simde__m128i mask = simde_mm_loadu_si128(
              reinterpret_cast<const simde__m128i *>(lut.masks.data() + mask_offset));
          combined_mask = simde_mm_or_si128(combined_mask, mask);
        }
      }

      visibility = simde_mm_andnot_si128(combined_mask, visibility);

      if (simde_mm_testz_si128(visibility, visibility)) {
        return 0.0;
      }
    }

    idx += 4;
  }

  while (idx < end) {
    const auto &nbr       = ctx.neighbors.items[idx];
    const std::size_t j   = nbr.index;
    const double vx       = cx - ctx.coords.x(j);
    const double vy       = cy - ctx.coords.y(j);
    const double vz       = cz - ctx.coords.z(j);
    const double v_mag_sq = vx * vx + vy * vy + vz * vz;

    if (v_mag_sq == 0.0) {
      throw std::invalid_argument("SASA invalid geometry: coincident atom centers.");
    }

    const double inv_dist      = fast_rsqrt_scalar(v_mag_sq);
    const double limit         = (nbr.threshold_sq - v_mag_sq - r_i_sq) * inv_two_r;
    const double cos_threshold = limit * inv_dist;

    if (cos_threshold >= 1.0) return 0.0;

    if (cos_threshold > -1.0) {
      const double nx               = vx * inv_dist;
      const double ny               = vy * inv_dist;
      const double nz               = vz * inv_dist;
      const int dir_idx             = dir_bin_from_unit(nx, ny, nz, lut.dir_resolution);
      const int angle_idx           = angle_bin_from_cos(cos_threshold, lut.angle_bins);
      const std::size_t mask_offset = (static_cast<std::size_t>(dir_idx) * lut.angle_bins + angle_idx) * 2;
      const simde__m128i mask       = simde_mm_loadu_si128(
          reinterpret_cast<const simde__m128i *>(lut.masks.data() + mask_offset));
      visibility = simde_mm_andnot_si128(mask, visibility);

      if (simde_mm_testz_si128(visibility, visibility)) {
        return 0.0;
      }
    }
    ++idx;
  }

  alignas(16) std::uint64_t final_vis[2];
  simde_mm_store_si128(reinterpret_cast<simde__m128i *>(final_vis), visibility);
  const std::size_t accessible = popcount64(final_vis[0]) + popcount64(final_vis[1]);

  return (4.0 * Pi * r_i_sq / 128.0) * static_cast<double>(accessible);
}

// Batched kernel for 256-point bitmask (4 words).
template <typename CoordAccessor>
[[nodiscard]] double compute_bitmask_256_batched(const BitmaskLut &lut,
                                                 const ComputeContext<CoordAccessor> &ctx) {
  const double r_i        = ctx.radius();
  const double r_i_sq     = ctx.radius_sq();
  const std::size_t start = ctx.neighbor_start();
  const std::size_t end   = ctx.neighbor_end();

  if (start == end) return 4.0 * Pi * r_i_sq;

  const double inv_two_r = 1.0 / (2.0 * r_i);
  const double cx        = ctx.coords.x(ctx.atom_index);
  const double cy        = ctx.coords.y(ctx.atom_index);
  const double cz        = ctx.coords.z(ctx.atom_index);

  alignas(32) std::uint64_t vis_arr[4] = {~std::uint64_t{0},
                                          ~std::uint64_t{0},
                                          ~std::uint64_t{0},
                                          lut.last_word_mask};
  simde__m256i visibility = simde_mm256_load_si256(reinterpret_cast<const simde__m256i *>(vis_arr));

  const simde__m256d cx_vec        = simde_mm256_set1_pd(cx);
  const simde__m256d cy_vec        = simde_mm256_set1_pd(cy);
  const simde__m256d cz_vec        = simde_mm256_set1_pd(cz);
  const simde__m256d r_i_sq_vec    = simde_mm256_set1_pd(r_i_sq);
  const simde__m256d inv_two_r_vec = simde_mm256_set1_pd(inv_two_r);
  const simde__m256d neg_one       = simde_mm256_set1_pd(-1.0);
  const simde__m256d pos_one       = simde_mm256_set1_pd(1.0);

  std::size_t idx = start;

  while (idx + 4 <= end) {
    alignas(32) double nbr_cx[4];
    alignas(32) double nbr_cy[4];
    alignas(32) double nbr_cz[4];
    alignas(32) double nbr_thresh_sq[4];

    for (int b = 0; b < 4; ++b) {
      const auto &nbr     = ctx.neighbors.items[idx + b];
      const std::size_t j = nbr.index;
      nbr_cx[b]           = ctx.coords.x(j);
      nbr_cy[b]           = ctx.coords.y(j);
      nbr_cz[b]           = ctx.coords.z(j);
      nbr_thresh_sq[b]    = nbr.threshold_sq;
    }

    simde__m256d nbr_cx_vec    = simde_mm256_load_pd(nbr_cx);
    simde__m256d nbr_cy_vec    = simde_mm256_load_pd(nbr_cy);
    simde__m256d nbr_cz_vec    = simde_mm256_load_pd(nbr_cz);
    simde__m256d thresh_sq_vec = simde_mm256_load_pd(nbr_thresh_sq);

    simde__m256d vx    = simde_mm256_sub_pd(cx_vec, nbr_cx_vec);
    simde__m256d vy    = simde_mm256_sub_pd(cy_vec, nbr_cy_vec);
    simde__m256d vz    = simde_mm256_sub_pd(cz_vec, nbr_cz_vec);
    simde__m256d vx_sq = simde_mm256_mul_pd(vx, vx);
    simde__m256d vy_sq = simde_mm256_mul_pd(vy, vy);
    simde__m256d vz_sq = simde_mm256_mul_pd(vz, vz);

    simde__m256d v_mag_sq  = simde_mm256_add_pd(vx_sq, simde_mm256_add_pd(vy_sq, vz_sq));
    simde__m256d numerator = simde_mm256_sub_pd(simde_mm256_sub_pd(thresh_sq_vec, v_mag_sq), r_i_sq_vec);
    simde__m256d limit     = simde_mm256_mul_pd(numerator, inv_two_r_vec);

    simde__m256d inv_dist      = fast_rsqrt_avx2(v_mag_sq);
    simde__m256d cos_threshold = simde_mm256_mul_pd(limit, inv_dist);

    simde__m256d full_occlude_mask = simde_mm256_cmp_pd(cos_threshold, pos_one, SIMDE_CMP_GE_OQ);
    if (simde_mm256_movemask_pd(full_occlude_mask) != 0) {
      alignas(32) double cos_arr[4];
      simde_mm256_store_pd(cos_arr, cos_threshold);
      for (int b = 0; b < 4; ++b) {
        if (cos_arr[b] >= 1.0) {
          return 0.0;
        }
      }
    }

    simde__m256d valid_mask = simde_mm256_cmp_pd(cos_threshold, neg_one, SIMDE_CMP_GT_OQ);
    int valid_bits          = simde_mm256_movemask_pd(valid_mask);

    if (valid_bits != 0) {
      simde__m256d nx = simde_mm256_mul_pd(vx, inv_dist);
      simde__m256d ny = simde_mm256_mul_pd(vy, inv_dist);
      simde__m256d nz = simde_mm256_mul_pd(vz, inv_dist);

      simde__m256d u, v;
      octa_encode_avx2(nx, ny, nz, u, v);

      alignas(16) int dir_indices[4];
      alignas(16) int angle_indices[4];
      dir_bin_avx2(u, v, lut.dir_resolution, dir_indices);
      angle_bin_avx2(cos_threshold, lut.angle_bins, angle_indices);

      simde__m256i combined_mask = simde_mm256_setzero_si256();
      for (int b = 0; b < 4; ++b) {
        if (valid_bits & (1 << b)) {
          const std::size_t mask_offset = (static_cast<std::size_t>(dir_indices[b]) * lut.angle_bins +
                                           angle_indices[b]) *
                                          4;
          const simde__m256i mask = simde_mm256_loadu_si256(
              reinterpret_cast<const simde__m256i *>(lut.masks.data() + mask_offset));
          combined_mask = simde_mm256_or_si256(combined_mask, mask);
        }
      }

      visibility = simde_mm256_andnot_si256(combined_mask, visibility);

      if (simde_mm256_testz_si256(visibility, visibility)) {
        return 0.0;
      }
    }

    idx += 4;
  }

  while (idx < end) {
    const auto &nbr       = ctx.neighbors.items[idx];
    const std::size_t j   = nbr.index;
    const double vx       = cx - ctx.coords.x(j);
    const double vy       = cy - ctx.coords.y(j);
    const double vz       = cz - ctx.coords.z(j);
    const double v_mag_sq = vx * vx + vy * vy + vz * vz;

    if (v_mag_sq == 0.0) {
      throw std::invalid_argument("SASA invalid geometry: coincident atom centers.");
    }

    const double inv_dist      = fast_rsqrt_scalar(v_mag_sq);
    const double limit         = (nbr.threshold_sq - v_mag_sq - r_i_sq) * inv_two_r;
    const double cos_threshold = limit * inv_dist;

    if (cos_threshold >= 1.0) return 0.0;

    if (cos_threshold > -1.0) {
      const double nx               = vx * inv_dist;
      const double ny               = vy * inv_dist;
      const double nz               = vz * inv_dist;
      const int dir_idx             = dir_bin_from_unit(nx, ny, nz, lut.dir_resolution);
      const int angle_idx           = angle_bin_from_cos(cos_threshold, lut.angle_bins);
      const std::size_t mask_offset = (static_cast<std::size_t>(dir_idx) * lut.angle_bins + angle_idx) * 4;
      const simde__m256i mask       = simde_mm256_loadu_si256(
          reinterpret_cast<const simde__m256i *>(lut.masks.data() + mask_offset));
      visibility = simde_mm256_andnot_si256(mask, visibility);

      if (simde_mm256_testz_si256(visibility, visibility)) {
        return 0.0;
      }
    }
    ++idx;
  }

  alignas(32) std::uint64_t final_vis[4];
  simde_mm256_store_si256(reinterpret_cast<simde__m256i *>(final_vis), visibility);
  const std::size_t accessible = popcount64(final_vis[0]) + popcount64(final_vis[1]) +
                                 popcount64(final_vis[2]) + popcount64(final_vis[3]);

  return (4.0 * Pi * r_i_sq / 256.0) * static_cast<double>(accessible);
}

} // namespace lahuta::analysis::simd

#endif // LAHUTA_ANALYSIS_SASA_METHODS_BITMASK_SIMD_HPP
