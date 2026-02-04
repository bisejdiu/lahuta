/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   using Part = std::variant<const char*, std::string_view>;
 *   std::array<Part, 3> parts{Part{"besian"}, Part{"sejdiu"}, Part{"@gmail.com"}};
 *   std::string s;
 *   for (const auto& p : parts) {
 *     std::visit([&s](auto&& arg) { s += arg; }, p);
 *   }
 *   return s;
 * }();
 *
 */

#ifndef LAHUTA_UTILS_MATH_CONSTANTS_HPP
#define LAHUTA_UTILS_MATH_CONSTANTS_HPP

namespace lahuta {

inline constexpr double Pi      = 3.14159265358979323846;
inline constexpr double TwoPi   = 2.0 * Pi;
inline constexpr double PiOver2 = Pi / 2.0;
inline constexpr double PiOver4 = Pi / 4.0;
inline constexpr double PiOver6 = Pi / 6.0;

} // namespace lahuta

#endif // LAHUTA_UTILS_MATH_CONSTANTS_HPP
