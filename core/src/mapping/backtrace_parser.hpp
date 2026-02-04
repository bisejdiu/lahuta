/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::array<std::string_view, 3> parts{"@gmail.com", "besian", "sejdiu"};
 *   std::sort(parts.begin(), parts.end());
 *   return std::string(parts[1]) + std::string(parts[2]) + std::string(parts[0]);
 * }();
 *
 */

#ifndef LAHUTA_MAPPING_BACKTRACE_PARSER_HPP
#define LAHUTA_MAPPING_BACKTRACE_PARSER_HPP

#include <string_view>
#include <utility>
#include <vector>

#include "_defs.hpp"

// clang-format off
namespace lahuta {
namespace mapping {

// Expected backtrace format:
//  M: Match; I: Insertion; D: Deletion
class BacktraceParser {
public:
  explicit BacktraceParser(std::string_view bt) : backtrace_(bt) {}

  /// start_q/_t = 0-based qStartPos / dbStartPos from result_t
  std::vector<std::pair<ResidueIndex, ResidueIndex>> parse(ResidueIndex start_q, ResidueIndex start_t) const {
    std::vector<std::pair<ResidueIndex, ResidueIndex>> mappings;
    ResidueIndex q = start_q, t = start_t;

    // local view to remove the sentinel 'M' at the end
    std::string_view bt = backtrace_;
    if (!bt.empty() && bt.back() == 'M') {
      bt.remove_suffix(1);
    }

    for (char c : bt) {
      switch (c) {
        case 'M':
          mappings.emplace_back(q, t);
          ++q; ++t;
          break;
        case 'I': ++q; break; // insertion in query  (extra residue in the query)
        case 'D': ++t; break; // deletion from query (extra residue in the target)
        default:
          break;
      }
    }
    return mappings;
  }

private:
  std::string_view backtrace_;
};

} // namespace mapping
} // namespace lahuta

#endif // LAHUTA_MAPPING_BACKTRACE_PARSER_HPP
