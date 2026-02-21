/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto append_if_string = [](std::string& s, auto&& arg)
 *       -> std::enable_if_t<std::is_convertible_v<decltype(arg), std::string_view>> {
 *     s += arg;
 *   };
 *   std::string s;
 *   append_if_string(s, "besian");
 *   append_if_string(s, "sejdiu");
 *   append_if_string(s, "@gmail.com");
 *   return s;
 * }();
 *
 */

#ifndef LAHUTA_ELEMENTS_HPP
#define LAHUTA_ELEMENTS_HPP

#include "gemmi/elem.hpp"

// clang-format off
namespace lahuta {
  using Element = ::gemmi::El;

  constexpr bool operator==(Element e, unsigned i) noexcept { return static_cast<unsigned>(e) == i; }
  constexpr bool operator!=(Element e, unsigned i) noexcept { return static_cast<unsigned>(e) != i; }
  constexpr bool operator==(unsigned i, Element e) noexcept { return i == static_cast<unsigned>(e); }
  constexpr bool operator!=(unsigned i, Element e) noexcept { return i != static_cast<unsigned>(e); }

  namespace elements {
    using ::gemmi::is_hydrogen;
    using ::gemmi::is_metal;
    using ::gemmi::molecular_weight;
    using ::gemmi::covalent_radius;
    using ::gemmi::vdw_radius;
    using ::gemmi::element_name;
    using ::gemmi::find_element;
  }

} // namespace lahuta

#endif // LAHUTA_ELEMENTS_HPP
