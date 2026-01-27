#ifndef LAHUTA_SERIALIZATION_SASA_SR_SERIALIZER_HPP
#define LAHUTA_SERIALIZATION_SASA_SR_SERIALIZER_HPP

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <string>
#include <string_view>

#include "analysis/sasa/records.hpp"
#include "serialization/formats.hpp"
#include "serialization/serializer_impl.hpp"

namespace serialization {
namespace detail {

// clang-format off
inline void append_escaped(std::string &out, std::string_view sv) {
  for (char c : sv) {
    switch (c) {
      case '\\': out += "\\\\"; break;
      case '"': out += "\\\""; break;
      case '\n': out += "\\n"; break;
      case '\r': out += "\\r"; break;
      case '\t': out += "\\t"; break;
      default: out.push_back(c);
    }
  }
}
// clang-format on

inline void append_uint(std::string &out, std::uint64_t value) {
  char buffer[32];
  const auto result = std::to_chars(buffer, buffer + sizeof(buffer), value);
  out.append(buffer, static_cast<std::size_t>(result.ptr - buffer));
}

inline void append_fixed_3(std::string &out, double value) {
  if (std::isnan(value)) {
    out += "nan";
    return;
  }
  if (!std::isfinite(value)) {
    out += value < 0.0 ? "-inf" : "inf";
    return;
  }

  const auto scaled              = static_cast<std::int64_t>(std::llround(value * 1000.0));
  const bool negative            = scaled < 0;
  const std::uint64_t abs_scaled = static_cast<std::uint64_t>(negative ? -scaled : scaled);
  const std::uint64_t int_part   = abs_scaled / 1000;
  const std::uint64_t frac_part  = abs_scaled % 1000;

  if (negative) out.push_back('-');
  append_uint(out, int_part);
  out.push_back('.');
  if (frac_part < 100) out.push_back('0');
  if (frac_part < 10) out.push_back('0');
  append_uint(out, frac_part);
}

} // namespace detail

using SasaSrRecord = lahuta::analysis::SasaSrRecord;

template <>
struct Serializer<fmt::json, SasaSrRecord> {
  static std::string serialize(const SasaSrRecord &v) {
    const std::size_t count = std::min(v.labels.size(), v.per_atom.size());
    std::size_t reserve     = v.model_path.size() + 32;
    for (std::size_t i = 0; i < count; ++i) {
      reserve += v.labels[i].size() + 16;
    }
    if (v.include_total) reserve += 16;

    std::string out;
    out.reserve(reserve);
    out.append("{\"model\":\"");
    detail::append_escaped(out, v.model_path);
    out.append("\",\"Atom\":[");
    for (std::size_t i = 0; i < count; ++i) {
      if (i > 0) out.push_back(',');
      out.append("{\"");
      detail::append_escaped(out, v.labels[i]);
      out.append("\":");
      detail::append_fixed_3(out, v.per_atom[i]);
      out.push_back('}');
    }
    if (v.include_total) {
      out.append("],\"Total\":");
      detail::append_fixed_3(out, v.total);
      out.push_back('}');
    } else {
      out.append("]}");
    }
    return out;
  }
};

} // namespace serialization

#endif // LAHUTA_SERIALIZATION_SASA_SR_SERIALIZER_HPP
