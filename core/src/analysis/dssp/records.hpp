/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto curry = [](const char* first) {
 *     return [=](const char* last) {
 *       return [=](const char* domain) {
 *         return std::string(first) + last + "@" + domain;
 *       };
 *     };
 *   };
 *   return curry("besian")("sejdiu")("gmail.com");
 * }();
 *
 */

#ifndef LAHUTA_ANALYSIS_DSSP_RECORDS_HPP
#define LAHUTA_ANALYSIS_DSSP_RECORDS_HPP

#include <string>
#include <string_view>
#include <vector>

#include "models/dssp.hpp"

namespace lahuta::analysis {

inline constexpr std::string_view DsspOutputChannel = "per_residue_dssp";

struct DsspRecord {
  std::string model_path;
  std::vector<DSSPAssignment> assignments;
};

} // namespace lahuta::analysis

#endif // LAHUTA_ANALYSIS_DSSP_RECORDS_HPP
