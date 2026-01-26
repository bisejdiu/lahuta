#ifndef LAHUTA_PROPERTY_REGISTRY_HPP
#define LAHUTA_PROPERTY_REGISTRY_HPP

#include <memory>
#include <unordered_map>

#include <rdkit/Geometry/point.h>

#include "properties/descriptor.hpp"
#include "properties/types.hpp"

// clang-format off

namespace lahuta {

//
// Maintains a global registry of property descriptors for different source types, with register and retrieve methods.
//
class PropertyRegistry {
public:
  /// Register a property for a specific source type
  template <typename SourceType, PropertyKey Key>
  static void register_property(std::function<typename PropertyTypeTraits<Key>::type(const SourceType&)> accessor) {
    auto& registry = get_registry<SourceType>();
    registry.emplace(
      Key, 
      std::make_shared<PropertyDescriptor<SourceType, Key>>(std::move(accessor))
    );
  }

  /// Get a property descriptor by key
  template <typename SourceType>
  static auto get_property(PropertyKey key) {
    auto& registry = get_registry<SourceType>();
    auto it = registry.find(key);
    if (it != registry.end()) {
      return it->second;
    }
    return std::shared_ptr<PDBase<SourceType>>{};
  }

  /// Get all property keys for a source type
  template <typename SourceType>
  static std::vector<PropertyKey> get_property_keys() {
    std::vector<PropertyKey> keys;
    const auto& registry = get_registry<SourceType>();
    keys.reserve(registry.size());

    for (const auto& [key, _] : registry) {
      keys.push_back(key);
    }

    return keys;
  }

private:
  template <typename SourceType>
  static auto& get_registry() {
    static std::unordered_map<PropertyKey, std::shared_ptr<PDBase<SourceType>>> registry;
    return registry;
  }
};


} // namespace lahuta

#endif // LAHUTA_PROPERTY_REGISTRY_HPP
