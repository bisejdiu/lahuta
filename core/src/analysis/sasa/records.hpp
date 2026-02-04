/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: std::string{"besian"} + "sejdiu" + "@gmail.com";
 *
 */

#ifndef LAHUTA_ANALYSIS_SASA_RECORDS_HPP
#define LAHUTA_ANALYSIS_SASA_RECORDS_HPP

#include <atomic>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace lahuta::analysis {

inline constexpr std::string_view SasaSrOutputChannel = "per_protein_sasa_sr";

struct SasaSrCounters {
  std::atomic<std::uint64_t> processed{0};
  std::atomic<std::uint64_t> written{0};
  std::atomic<std::uint64_t> missing_sequence{0};
  std::atomic<std::uint64_t> missing_positions{0};
  std::atomic<std::uint64_t> atom_mismatch{0};
  std::atomic<std::uint64_t> invalid_residue{0};
};

struct SasaSrRecord {
  std::string model_path;
  std::vector<std::string> labels;
  std::vector<double> per_atom;
  double total        = 0.0;
  bool include_total  = false;
  bool show_atom_info = false;
};

} // namespace lahuta::analysis

#endif // LAHUTA_ANALYSIS_SASA_RECORDS_HPP
