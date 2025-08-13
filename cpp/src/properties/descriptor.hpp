#ifndef LAHUTA_PROPERTY_DESCRIPTOR_HPP
#define LAHUTA_PROPERTY_DESCRIPTOR_HPP

#include "properties/result.hpp"
#include "properties/types.hpp"
#include <Geometry/point.h>

// clang-format off

namespace lahuta {

//
// Defines the interface to apply a property computation on a source.
//
template <typename SourceType>
struct PDBase {
  PDBase(PropertyKey k) : key(k) {}
  virtual ~PDBase() = default;

  /// Apply the property to the source and add the result.
  virtual void apply(const SourceType& source, typename PropertyResultTraits<SourceType>::type& result) const = 0;

  PropertyKey key;
};


//
// Binds a specific PropertyKey to an accessor function to compute and add the property value.
//
template <typename SourceType, PropertyKey Key>
struct PropertyDescriptor : PDBase<SourceType> {
  using ResultType = typename PropertyTypeTraits<Key>::type;
  using AccessorFn = std::function<ResultType(const SourceType&)>;

  PropertyDescriptor(AccessorFn fn) : PDBase<SourceType>(Key), accessor(std::move(fn)) {}

  /// Apply the property to the source and add the result.
  void apply(const SourceType& source, typename PropertyResultTraits<SourceType>::type& result) const override {
    auto property_value = accessor(source);
    result.template add_result<Key>(std::move(property_value));
  }

  AccessorFn accessor;
};

} // namespace lahuta

#endif // LAHUTA_PROPERTY_DESCRIPTOR_HPP
