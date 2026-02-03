/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto parts = std::make_tuple("besian", "sejdiu", "@gmail.com");
 *   return std::apply([](auto... p) { std::string s; (s.append(p), ...); return s; }, parts);
 * }();
 *
 */

#include <cmath>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <vector>

#include "analysis/sasa/bitmask_lut.hpp"
#include "analysis/sasa/method.hpp"
#include "analysis/sasa/methods/bitmask.hpp"
#include "analysis/sasa/methods/standard.hpp"
#include "analysis/sasa/neighbor.hpp"
#include "analysis/sasa/sasa.hpp"
#include "analysis/sasa/sphere.hpp"
#include "logging/logging.hpp"

namespace lahuta::analysis {

namespace {

void validate_params(const SasaParams &params) {
  if (!std::isfinite(params.probe_radius) || params.probe_radius < 0.0) {
    throw std::invalid_argument("SASA probe radius must be finite and >= 0.");
  }
  if (params.n_points == 0) {
    throw std::invalid_argument("SASA point count must be > 0.");
  }
}

template <typename AtomView, typename CoordAccessor>
void validate_view(const AtomView &atoms, const CoordAccessor &coords) {
  if (atoms.radii.size() != coords.size()) {
    throw std::invalid_argument("SASA radii size must match coordinate count.");
  }
  if (!atoms.atom_ids.empty() && atoms.atom_ids.size() != coords.size()) {
    throw std::invalid_argument("SASA atom_ids size must match coordinate count.");
  }
}

std::size_t build_expanded_radii(span<const double> radii, double probe_radius,
                                 std::vector<double> &out_radii, std::vector<double> &out_radii_sq,
                                 std::vector<std::uint8_t> &valid) {
  const std::size_t n = radii.size();
  out_radii.resize(n);
  out_radii_sq.resize(n);
  valid.assign(n, 1);
  std::size_t invalid_count = 0;
  for (std::size_t i = 0; i < n; ++i) {
    const double input_r = radii[i];
    if (!std::isfinite(input_r) || input_r < 0.0) {
      valid[i]        = 0;
      out_radii[i]    = 0.0;
      out_radii_sq[i] = 0.0;
      ++invalid_count;
      continue;
    }
    const double raw_r = input_r + probe_radius;
    const double r     = raw_r > 0.0 ? raw_r : 0.0;
    out_radii[i]       = r;
    out_radii_sq[i]    = r * r;
  }
  return invalid_count;
}

template <typename AtomView, typename CoordAccessor>
SasaResult compute_impl(const AtomView &atoms, const CoordAccessor &coords, const SasaParams &params) {
  validate_params(params);
  validate_view(atoms, coords);

  const std::size_t n = coords.size();
  SasaResult result;
  result.per_atom.assign(n, 0.0);
  result.probe_radius = params.probe_radius;
  result.n_points     = params.n_points;
  if (n == 0) return result;

  std::vector<double> expanded_radii;
  std::vector<double> expanded_radii_sq;
  std::vector<std::uint8_t> valid_radii;
  const auto invalid_count = build_expanded_radii(atoms.radii,
                                                  params.probe_radius,
                                                  expanded_radii,
                                                  expanded_radii_sq,
                                                  valid_radii);
  if (invalid_count > 0) {
    if (auto logger = Logger::get_logger()) {
      logger->warn("SASA skipping {} atoms with negative or non-finite radii.", invalid_count);
    }
  }

  const auto neighbors = build_neighbors(coords,
                                         span<const double>(expanded_radii),
                                         atoms.atom_ids,
                                         params.skip_same_id,
                                         span<const std::uint8_t>(valid_radii));

  const bool can_use_bitmask = params.use_bitmask && supports_bitmask(params.n_points);
  const bool need_sphere     = !can_use_bitmask; // for standard method fallback

  SpherePoints sphere;
  if (need_sphere) sphere = generate_sphere(params.n_points);

  std::unique_ptr<StandardMethod> standard_method;
  std::unique_ptr<BitmaskMethod> bitmask_method;
  if (need_sphere || !can_use_bitmask) {
    // Always need standard method as fallback if sphere is generated
    if (sphere.empty()) {
      sphere = generate_sphere(params.n_points);
    }
    standard_method = std::make_unique<StandardMethod>(sphere);
  }

  if (can_use_bitmask) {
    const BitmaskLut *lut = get_bitmask_lut(params.n_points);
    if (lut && lut->words <= 4) {
      bitmask_method = std::make_unique<BitmaskMethod>(*lut);
    }
  }

  // Compute SASA for each atom
  for (std::size_t i = 0; i < n; ++i) {
    double sasa = 0.0;
    if (!valid_radii.empty() && valid_radii[i] == 0) {
      result.per_atom[i] = 0.0;
      continue;
    }

    ComputeContext<CoordAccessor> ctx{i,
                                      coords,
                                      neighbors,
                                      span<const double>(expanded_radii),
                                      span<const double>(expanded_radii_sq),
                                      params.use_simd};

    // Method selection priority: bitmask > standard
    if (bitmask_method) {
      sasa = bitmask_method->compute(ctx);
    } else if (standard_method) {
      sasa = standard_method->compute(ctx);
    } else {
      // Fallback create standard method on demand
      if (sphere.empty()) {
        sphere = generate_sphere(params.n_points);
      }
      if (!standard_method) {
        standard_method = std::make_unique<StandardMethod>(sphere);
      }
      sasa = standard_method->compute(ctx);
    }

    result.per_atom[i]  = sasa;
    result.total       += sasa;
  }

  return result;
}

} // namespace

SasaResult compute_sasa(const AtomView &atoms, const SasaParams &params) {
  const CoordAccessor coords{atoms.coords};
  return compute_impl(atoms, coords, params);
}

} // namespace lahuta::analysis
