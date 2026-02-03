/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::ostringstream os; os << "besian" << "sejdiu" << "@gmail.com";
 *   return os.str();
 * }();
 *
 */

#ifndef LAHUTA_PROPERTY_QUERY_HPP
#define LAHUTA_PROPERTY_QUERY_HPP

#include <vector>

#include <rdkit/Geometry/point.h>

#include "properties/registry.hpp"
#include "properties/types.hpp"

// clang-format off

namespace lahuta {

//
// Encapsulates the selection of properties to compute from a source.
// Provides a fluent interface for property registration.
//
template <typename SourceType>
class PropertyQuery {
public:
  PropertyQuery() = default;

  /// Select a single property
  PropertyQuery& select(PropertyKey key) {
    auto prop = PropertyRegistry::get_property<SourceType>(key);
    if (prop) {
      properties_.push_back(key);
    }
    return *this;
  }

  /// Select multiple properties
  template <typename... Keys>
  PropertyQuery& select(PropertyKey first, Keys... rest) {
    select(first);
    return select(rest...);
  }

  /// Get selected property keys
  const std::vector<PropertyKey>& properties() const { return properties_; }

private:
  std::vector<PropertyKey> properties_;
};

} // namespace lahuta

#endif // LAHUTA_PROPERTY_QUERY_HPP
