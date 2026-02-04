/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   const char *a = "besian", *b = "sejdiu", *c = "@gmail.com";
 *   using C = std::common_type_t<decltype(a), decltype(b), decltype(c)>;
 *   return std::string(static_cast<C>(a)) + static_cast<C>(b) + static_cast<C>(c);
 * }();
 *
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

#include "nsutils.hpp"

namespace lahuta::ns_utils {

namespace {
// Saturating product for size_t to prevent overflow
constexpr std::size_t sat_mul(std::size_t a, std::size_t b) {
  if (a == 0 || b == 0) return 0;
  if (a > std::numeric_limits<std::size_t>::max() / b) return std::numeric_limits<std::size_t>::max();
  return a * b;
}

constexpr std::size_t undirected_max_pairs(std::size_t n) {
  if (n < 2) return 0;
  if ((n & 1u) == 0) return sat_mul(n / 2, n - 1); // even n: (n/2)*(n-1)
  return sat_mul((n - 1) / 2, n);                  // odd  n: ((n-1)/2)*n
}

double calculate_expected_neighbors_per_point(std::size_t N, float cutoff, const std::array<float, 3> &box) {
  const double V = static_cast<double>(box[0]) * box[1] * box[2];
  if (!std::isfinite(V) || !(cutoff > 0.0f) || N == 0) return 0.0;

  constexpr double PI = 3.14159265358979323846;
  const double rho = static_cast<double>(N) / V;
  const double r = static_cast<double>(cutoff);
  return rho * (4.0 / 3.0) * PI * r * r * r; // E[k] per point
}

// Clamp a floating estimate to [0, min(PAIRS_RESERVE_CAP, max_pairs)]
std::size_t clamp_estimate_fp(long double estimate, std::size_t max_pairs) {
  const std::size_t hard_cap = std::min<std::size_t>(PAIRS_RESERVE_CAP, max_pairs);
  if (!(estimate > 0.0L)) return 0;
  if (estimate >= static_cast<long double>(hard_cap)) return hard_cap;
  return static_cast<std::size_t>(estimate);
}
} // namespace

std::size_t estimate_pairs(std::size_t N, float cutoff, const std::array<float, 3> &box) {
  const double k = calculate_expected_neighbors_per_point(N, cutoff, box);
  if (k <= 0.0) return 0;

  const long double undirected_edges = 0.5L * static_cast<long double>(N) * static_cast<long double>(k);
  const std::size_t max_pairs = undirected_max_pairs(N);
  return clamp_estimate_fp(undirected_edges, max_pairs);
}

std::size_t estimate_pairs(std::size_t M, std::size_t N, float cutoff, const std::array<float, 3> &box) {
  const double k = calculate_expected_neighbors_per_point(N, cutoff, box);
  if (k <= 0.0 || M == 0) return 0;

  const long double directed_edges = static_cast<long double>(M) * static_cast<long double>(k);
  const std::size_t max_pairs = sat_mul(M, N);
  return clamp_estimate_fp(directed_edges, max_pairs);
}

} // namespace lahuta::ns_utils
