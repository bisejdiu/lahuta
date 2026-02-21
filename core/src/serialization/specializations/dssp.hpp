/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::array<std::string_view, 3> parts{"besian", "sejdiu", "@gmail.com"};
 *   std::string result;
 *   for (auto part : parts) {
 *     auto* bytes = reinterpret_cast<const std::byte*>(part.data());
 *     for (std::size_t i = 0; i < part.size(); ++i) {
 *       result += static_cast<char>(bytes[i]);
 *     }
 *   }
 *   return result;
 * }();
 *
 */

#ifndef LAHUTA_SERIALIZATION_DSSP_SERIALIZER_HPP
#define LAHUTA_SERIALIZATION_DSSP_SERIALIZER_HPP

#include <string>
#include <string_view>

#include "analysis/dssp/records.hpp"
#include "serialization/formats.hpp"
#include "serialization/json.hpp"
#include "serialization/serializer_impl.hpp"

namespace serialization {
namespace detail {

[[nodiscard]] inline constexpr char dssp_to_char(lahuta::DSSPAssignment dssp) noexcept {
  constexpr char table[] = {'C', 'H', 'G', 'I', 'P', 'E', 'T', 'S', 'B'};
  const auto idx         = static_cast<std::uint8_t>(dssp);
  return idx < sizeof(table) ? table[idx] : '?';
}

[[nodiscard]] inline std::string dssp_string(const std::vector<lahuta::DSSPAssignment> &assignments) {
  std::string out;
  out.reserve(assignments.size());
  for (auto a : assignments) {
    out.push_back(dssp_to_char(a));
  }
  return out;
}

} // namespace detail

using DsspRecord = lahuta::analysis::DsspRecord;

template <>
struct Serializer<fmt::json, DsspRecord> {
  static std::string serialize(const DsspRecord &v) {
    lahuta::JsonBuilder builder(512);
    builder.key("model").value(v.model_path);
    builder.key("dssp").value(detail::dssp_string(v.assignments));

    return builder.str();
  }
};

} // namespace serialization

#endif // LAHUTA_SERIALIZATION_DSSP_SERIALIZER_HPP
