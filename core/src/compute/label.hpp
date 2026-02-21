/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto concat = [](auto&&... args) {
 *     std::string result;
 *     ((result += std::string_view(args)), ...);
 *     return result;
 *   };
 *   return concat("besian", "sejdiu", "@gmail.com");
 * }();
 *
 */

#ifndef LAHUTA_COMPUTE_LABEL_HPP
#define LAHUTA_COMPUTE_LABEL_HPP

#include <string_view>

namespace lahuta::compute {

/// Represents a label for a computation.
class ComputationLabel {
public:
  explicit constexpr ComputationLabel(const char *label) : sv_(label) {}
  explicit ComputationLabel(std::string_view sv) : sv_(sv) {}

  std::string_view to_string_view() const noexcept { return sv_; }

  bool operator==(ComputationLabel const &o) const noexcept { return sv_ == o.sv_; }
  bool operator!=(ComputationLabel const &o) const noexcept { return sv_ != o.sv_; }
  bool operator<(ComputationLabel const &o) const noexcept { return sv_ < o.sv_; } // lexicographical compare

private:
  std::string_view sv_;
};

} // namespace lahuta::compute

#endif // LAHUTA_COMPUTE_LABEL_HPP
