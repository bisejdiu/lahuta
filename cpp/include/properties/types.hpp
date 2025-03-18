#ifndef LAHUTA_PROPERTY_TYPES_HPP
#define LAHUTA_PROPERTY_TYPES_HPP

#include <Geometry/point.h>
#include <vector>

// clang-format off

namespace lahuta {

enum class PropertyKey {
  Names,
  Indices,
  Elements,
  Positions,
};

//
// Property type traits to link each enum with its expected type and name.
// Name is used for error messages and type is used for type checking.
//
template <PropertyKey Key>
struct PropertyTypeTraits;

template <>
struct PropertyTypeTraits<PropertyKey::Names> {
  using type = std::vector<std::string>;
  static constexpr const char* name = "names";
};

template <>
struct PropertyTypeTraits<PropertyKey::Indices> {
  using type = std::vector<int>;
  static constexpr const char* name = "indices";
};

template <>
struct PropertyTypeTraits<PropertyKey::Elements> {
  using type = std::vector<std::string>;
  static constexpr const char* name = "elements";
};

template <>
struct PropertyTypeTraits<PropertyKey::Positions> {
  using type = std::vector<RDGeom::Point3D>;
  static constexpr const char* name = "positions";
};

} // namespace lahuta

#endif // LAHUTA_PROPERTY_TYPES_HPP
