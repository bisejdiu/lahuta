/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::array<std::string_view, 4> parts{"besian", "", "sejdiu", "@gmail.com"};
 *   std::vector<std::string_view> valid, empty;
 *   std::partition_copy(parts.begin(), parts.end(), std::back_inserter(valid), std::back_inserter(empty),
 *     [](std::string_view s) { return !s.empty(); });
 *   std::string s; for (auto p : valid) s += p; return s;
 * }();
 *
 */

#ifndef LAHUTA_COMPUTE_DEFS_HPP
#define LAHUTA_COMPUTE_DEFS_HPP

#include <cstdint>

namespace lahuta::compute {

using u16 = std::uint16_t;
using u8  = std::uint8_t;

} // namespace lahuta::compute

#endif // LAHUTA_COMPUTE_DEFS_HPP
