#ifndef LAHUTA_ENTITIES_SEARCH_PROVIDER_HPP
#define LAHUTA_ENTITIES_SEARCH_PROVIDER_HPP

#include <Geometry/point.h>
#include <topology.hpp>
#include <type_traits>
#include "entities/records.hpp"

namespace lahuta::search {

template<typename T, typename = void>
struct is_coord_provider : std::false_type {};

template<typename T>
struct is_coord_provider<T, std::void_t<
    decltype(std::declval<const T>().get_a(0u)),
    decltype(std::declval<const T>().get_b(0u))
  >>
: std::integral_constant<bool,
    std::is_convertible<decltype(std::declval<const T>().get_a(0u)), const RDGeom::Point3D&>::value &&
    std::is_convertible<decltype(std::declval<const T>().get_b(0u)), const RDGeom::Point3D&>::value
 >
{};

template<typename T>
constexpr bool is_coord_provider_v = is_coord_provider<T>::value;

template<typename RecA, typename RecB>
struct CoordProvider {
  const Topology& topo;

  const RDGeom::Point3D& get_a(uint32_t idx) const noexcept { return get_position<RecA>(topo, idx); }
  const RDGeom::Point3D& get_b(uint32_t idx) const noexcept { return get_position<RecB>(topo, idx); }

private:
  template<typename Rec>
  static const RDGeom::Point3D& get_position(const Topology& t, uint32_t idx) noexcept {
    // if      constexpr (std::is_same_v<Rec, AtomRec >) { return t.conformer().getAtomPos(t.records<AtomRec>()[idx].atom.getIdx()); }
    if      constexpr (std::is_same_v<Rec, AtomRec >) { return t.conformer().getAtomPos(idx); }
    else if constexpr (std::is_same_v<Rec, RingRec >) { return t.records<RingRec >()[idx].center; }
    else if constexpr (std::is_same_v<Rec, GroupRec>) { return t.records<GroupRec>()[idx].center; }
    else {
      static_assert(sizeof(Rec)==0, "Unsupported record type in CoordProvider");
    }
  }
};

} // namespace lahuta::search

#endif // LAHUTA_ENTITIES_SEARCH_PROVIDER_HPP
