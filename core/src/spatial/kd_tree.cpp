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

#include "kd_tree.hpp"

namespace lahuta {

template class KDTree3<float>;
template class KDTree3<double>;

} // namespace lahuta
