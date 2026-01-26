#ifndef LAHUTA_ENTITIES_SEARCH_PROVIDER_HPP
#define LAHUTA_ENTITIES_SEARCH_PROVIDER_HPP

#include <type_traits>

#include <rdkit/Geometry/point.h>

#include "compute/topology_snapshot.hpp"
#include "entities/records.hpp"
#include "topology.hpp"

namespace lahuta::search {
namespace C = lahuta::compute;

template <typename T, typename = void>
struct is_coord_provider : std::false_type {};

// clang-format off
template <typename T>
struct is_coord_provider<T, std::void_t<
    decltype(std::declval<const T>().get_a(0u)),
    decltype(std::declval<const T>().get_b(0u))
  >>
: std::integral_constant<bool,
    std::is_convertible<decltype(std::declval<const T>().get_a(0u)), const RDGeom::Point3D &>::value &&
    std::is_convertible<decltype(std::declval<const T>().get_b(0u)), const RDGeom::Point3D &>::value
 >
{};

// clang-format on
template <typename T>
constexpr bool is_coord_provider_v = is_coord_provider<T>::value;

template <typename RecA, typename RecB>
struct CoordProvider {
  const C::TopologySnapshot &ts;

  RDGeom::Point3D get_a(uint32_t idx) const noexcept { return get_position<RecA>(idx); }
  RDGeom::Point3D get_b(uint32_t idx) const noexcept { return get_position<RecB>(idx); }

private:
  template <typename Rec>
  RDGeom::Point3D get_position(uint32_t idx) const noexcept {
    if constexpr (std::is_same_v<Rec, AtomRec>) {
      return ts.conf.getAtomPos(idx);
    } else if constexpr (std::is_same_v<Rec, RingRec>) {
      return ts.topo.records<RingRec>()[idx].center(ts.conf);
    } else if constexpr (std::is_same_v<Rec, GroupRec>) {
      return ts.topo.records<GroupRec>()[idx].center(ts.conf);
    } else {
      static_assert(sizeof(Rec) == 0, "Unsupported record type in CoordProvider");
    }
  }
};

} // namespace lahuta::search

#endif // LAHUTA_ENTITIES_SEARCH_PROVIDER_HPP
