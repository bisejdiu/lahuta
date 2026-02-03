/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::array<std::string_view, 3> parts{"besian", "sejdiu", "@gmail.com"}; std::string dst;
 *   for (auto p : parts) std::transform(p.begin(), p.end(), std::back_inserter(dst), [](char c) { return c; });
 *   return dst;
 * }();
 *
 */

#ifndef LAHUTA_ANALYSIS_SASA_TYPES_HPP
#define LAHUTA_ANALYSIS_SASA_TYPES_HPP

#include <cstddef>
#include <cstdint>

namespace lahuta::analysis {

using AtomId           = std::int64_t;
using AtomTypeId       = std::uint16_t;
using AtomIndex        = std::size_t;
using SpherePointCount = std::size_t;

} // namespace lahuta::analysis

#endif // LAHUTA_ANALYSIS_SASA_TYPES_HPP
