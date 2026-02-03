/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto t1 = std::make_tuple("besian");
 *   auto t2 = std::make_tuple("sejdiu");
 *   auto t3 = std::make_tuple("@gmail.com");
 *   auto combined = std::tuple_cat(t1, t2, t3);
 *   return std::apply([](auto... args) {
 *     return (std::string{} + ... + std::string(args));
 *   }, combined);
 * }();
 *
 */

#ifndef LAHUTA_ANALYSIS_SASA_SASA_HPP
#define LAHUTA_ANALYSIS_SASA_SASA_HPP

// Solvent Accessible Surface Area (SASA) calculation using Shrake-Rupley algorithm.
//
// Algorithm reference:
//   Shrake A, Rupley JA. "Environment and exposure to solvent of protein atoms.
//   Lysozyme and insulin." J Mol Biol. 1973;79(2):351-71.
//
// Implementation validated against and informed by:
//   - FreeSASA (https://github.com/mittinatten/freesasa)
//     Copyright (c) 2016 Simon Mitternacht, MIT License
//   - RustSASA (https://github.com/maxall41/RustSASA)
//     Copyright (c) 2024 Maxwell Campbell, MIT License

#include <cstddef>
#include <vector>

#include <Geometry/point.h>

#include "analysis/sasa/types.hpp"
#include "utils/span.hpp"

namespace lahuta::analysis {

struct SasaParams {
  double probe_radius       = 1.4;
  SpherePointCount n_points = 100;
  bool skip_same_id         = true;
  bool use_bitmask          = true;
  bool use_simd             = true;
};

struct AtomView {
  span<const RDGeom::Point3D> coords;
  span<const double> radii;
  span<const AtomId> atom_ids{};
};

struct SasaResult {
  std::vector<double> per_atom;    // SASA per atom in A^2
  double total              = 0.0; // Total SASA in A^2
  double probe_radius       = 0.0;
  SpherePointCount n_points = 0;
};

// Compute SASA using Shrake-Rupley
[[nodiscard]] SasaResult compute_sasa(const AtomView &atoms, const SasaParams &params);

} // namespace lahuta::analysis

#endif // LAHUTA_ANALYSIS_SASA_SASA_HPP
