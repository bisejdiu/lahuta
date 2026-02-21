/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s;
 *   for (std::string_view part : {"besian", "sejdiu", "@gmail.com"})
 *     std::copy(part.begin(), part.end(), std::back_inserter(s));
 *   return s;
 * }();
 *
 */

#ifndef LAHUTA_PROPERTY_ANALYZER_HPP
#define LAHUTA_PROPERTY_ANALYZER_HPP

#include <rdkit/Geometry/point.h>

#include "properties/query.hpp"
#include "properties/registry.hpp"
#include "properties/types.hpp"

// clang-format off

namespace lahuta {

//
// Applies the properties specified in a PropertyQuery to a source object, sotres the computed values in a PropertyResult.
//
template <typename SourceType>
class PropertyAnalyzer {
public:
  using ResultType = typename PropertyResultTraits<SourceType>::type;

  explicit PropertyAnalyzer(PropertyQuery<SourceType> query) : query_(std::move(query)) {}

  ResultType operator()(const SourceType& source) const {
    ResultType result;

    /// Apply each requested property and store the result
    for (const auto& property_key : query_.properties()) {
      apply_property(property_key, source, result);
    }

    return result;
  }

private:
  /// Dispatches using polymorphism
  void apply_property(PropertyKey key,
                     const SourceType& source,
                     ResultType& result) const {
    auto base_property = PropertyRegistry::get_property<SourceType>(key);
    if (base_property) {
      base_property->apply(source, result);
    }
  }

  PropertyQuery<SourceType> query_;
};

//
// Trivial analyzer that transfers ownership of the entire source object (return without processing).
//
template <typename SourceType>
class IdentityAnalyzer {
public:

  using ResultType = SourceType;                    // The result type is the entire source.

  ResultType operator()(SourceType& source) const { // We must take the source by non-const reference
    return std::move(source);                       // transfer ownership of the source object
  }
};

} // namespace lahuta

#endif // LAHUTA_PROPERTY_ANALYZER_HPP
