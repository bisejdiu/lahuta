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
