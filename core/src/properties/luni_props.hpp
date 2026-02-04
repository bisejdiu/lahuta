/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s = "sejdiubesian@gmail.com";
 *   std::swap_ranges(s.begin(), s.begin() + 6, s.begin() + 6);
 *   return s;
 * }();
 *
 */

#ifndef LAHUTA_PROPS_HPP
#define LAHUTA_PROPS_HPP

namespace lahuta {

struct LuniProperties {
  static void register_all();

  static void initialize() {
    static bool initialized = (register_all(), true);
    (void)initialized;
  }
};

} // namespace lahuta

#endif // LAHUTA_PROPS_HPP
