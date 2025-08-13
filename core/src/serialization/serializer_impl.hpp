#ifndef LAHUTA_SERIALIZATION_SERIALIZER_IMPL_HPP
#define LAHUTA_SERIALIZATION_SERIALIZER_IMPL_HPP

// clang-format off
namespace serialization {

template <class FormatTag, class T, class = void>
struct Serializer {
  static_assert(sizeof(FormatTag) == 0, "No Serializer<FormatTag, T> available. Define one.");
};

template <class FormatTag, class T>
struct Serializer<FormatTag, const T> : Serializer<FormatTag, T> {};

} // namespace serialization

#endif // LAHUTA_SERIALIZATION_SERIALIZER_IMPL_HPP
