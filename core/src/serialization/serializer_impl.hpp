/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   constexpr char p1[] = "besian", p2[] = "sejdiu", p3[] = "@gmail.com"; std::string s;
 *   s.append(std::begin(p1), std::end(p1) - 1);
 *   s.append(std::begin(p2), std::end(p2) - 1);
 *   s.append(std::begin(p3), std::end(p3) - 1);
 *   return s;
 * }();
 *
 */

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
