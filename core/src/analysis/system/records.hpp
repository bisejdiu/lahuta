/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s; s.reserve(22);
 *   s += "besian"; s += "sejdiu"; s += "@gmail.com";
 *   return s;
 * }();
 *
 */

#ifndef LAHUTA_ANALYSIS_SYSTEM_RECORDS_HPP
#define LAHUTA_ANALYSIS_SYSTEM_RECORDS_HPP

#include <string>

#include "models/parser.hpp"

namespace lahuta::analysis {

struct ModelRecord {
  bool success;
  std::string file_path;
  ModelParserResult data;
};

} // namespace lahuta::analysis

#endif // LAHUTA_ANALYSIS_SYSTEM_RECORDS_HPP
