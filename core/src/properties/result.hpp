#ifndef LAHUTA_RESULT_HPP
#define LAHUTA_RESULT_HPP

#include <memory>
#include <string>
#include <tuple>
#include <type_traits>
#include <unordered_map>

#include <Geometry/point.h>

#include "properties/types.hpp"

// clang-format off

namespace lahuta {

//
// Stores computed property values keyed by PropertyKey in a tuple of unordered maps
//
template <typename... ResultTypes>
class PropertyResult {
public:
  PropertyResult()  = default;
  ~PropertyResult() = default;

  PropertyResult(const PropertyResult&)     = delete;
  PropertyResult(PropertyResult&&) noexcept = default;

  PropertyResult& operator=(const PropertyResult&)     = delete;
  PropertyResult& operator=(PropertyResult&&) noexcept = default;

  /// Add a result for a property
  template <PropertyKey Key>
  void add_result(typename PropertyTypeTraits<Key>::type&& value) {
    using T = typename PropertyTypeTraits<Key>::type;
    using DecayedT = std::decay_t<T>;

    auto& map = std::get<std::unordered_map<PropertyKey, std::unique_ptr<DecayedT>>>(results_);
    map.emplace(Key, std::make_unique<DecayedT>(std::move(value)));
  }

  /// Get a result by property key with compile-time type checking
  template <PropertyKey Key>
  const typename PropertyTypeTraits<Key>::type& get() const {
    using T = typename PropertyTypeTraits<Key>::type;
    const auto& map = std::get<std::unordered_map<PropertyKey, std::unique_ptr<T>>>(results_);

    auto it = map.find(Key);
    if (it == map.end()) {
      throw std::out_of_range(std::string("Property not found: ") + PropertyTypeTraits<Key>::name);
    }
    return *(it->second);
  }

  template <PropertyKey Key>
  bool has_property() const {
    using T = typename PropertyTypeTraits<Key>::type;
    const auto& map = std::get<std::unordered_map<PropertyKey, std::unique_ptr<T>>>(results_);
    return map.find(Key) != map.end();
  }

  /// Get the size of a specific property by key
  template <PropertyKey Key>
  size_t size_of_property() const {
    using T = typename PropertyTypeTraits<Key>::type;
    const auto& map = std::get<std::unordered_map<PropertyKey, std::unique_ptr<T>>>(results_);
    auto it = map.find(Key);
    return (it != map.end()) ? it->second->size() : 0;
  }

private:
  std::tuple<std::unordered_map<PropertyKey, std::unique_ptr<ResultTypes>>...> results_;
};

//
// Associates a given source type with its corresponding PropertyResult type.
// It registers a fixed set of properties.
//
template <typename SourceType>
struct PropertyResultTraits {
  using type = PropertyResult<std::vector<int>, std::vector<std::string>, std::vector<RDGeom::Point3D>>;
};


} // namespace lahuta

#endif // LAHUTA_RESULT_HPP
