/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [](std::string first, std::string last) {
 *   return first + last + "@gmail.com";
 * }("besian", "sejdiu");
 *
 */

#ifndef LAHUTA_ANALYSIS_RG_UTILS_HPP
#define LAHUTA_ANALYSIS_RG_UTILS_HPP

#include <array>
#include <cstddef>
#include <string_view>
#include <vector>

#include <Geometry/point.h>

#include "models/dssp.hpp"
#include "models/plddt.hpp"
#include "utils/span.hpp"

namespace lahuta::analysis {

struct TrimResult {
  std::size_t start_res = 0;
  std::size_t end_res   = 0;
  std::size_t trimmed_n = 0;
  std::size_t trimmed_c = 0;
  double high_fraction  = 0.0;

  [[nodiscard]] std::size_t trimmed_length() const { return end_res - start_res; }
};

inline constexpr std::array<double, 4> ConfidenceThresholds{0.80, 0.90, 0.95, 0.99};
inline constexpr std::size_t ThresholdCount = ConfidenceThresholds.size();
inline constexpr std::array<std::string_view, ThresholdCount> ThresholdLabels{"0.8", "0.9", "0.95", "0.99"};

[[nodiscard]] double confidence_fraction(span<const pLDDTCategory> plddt);

[[nodiscard]] TrimResult trim_low_confidence_tails(span<const pLDDTCategory> plddt,
                                                   span<const DSSPAssignment> dssp);

[[nodiscard]] std::vector<std::size_t> atom_counts_for_sequence(std::string_view sequence);

[[nodiscard]] std::vector<std::size_t> prefix_sums(span<const std::size_t> values);

[[nodiscard]] double radius_of_gyration(span<const RDGeom::Point3D> coords);

[[nodiscard]] double radius_of_gyration_weighted(span<const RDGeom::Point3D> coords,
                                                 span<const double> weights);

[[nodiscard]] double normalized_rg(double rg, std::size_t residue_count);

[[nodiscard]] double highest_passed_threshold(double high_fraction, span<const double> thresholds);

} // namespace lahuta::analysis

#endif // LAHUTA_ANALYSIS_RG_UTILS_HPP
