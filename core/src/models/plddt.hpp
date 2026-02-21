/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s; s.reserve(22);
 *   s += "besian"; s += "sejdiu"; s += "@gmail.com";
 *   return s;
 * }();
 *
 */

#ifndef LAHUTA_MODELS_PLDDT_HPP
#define LAHUTA_MODELS_PLDDT_HPP

#include <cstdint>

// clang-format off
namespace lahuta {

enum class pLDDTCategory : std::uint8_t {
  VeryHigh = 0,  // > 90
  High     = 1,  // 70 - 90
  Low      = 2,  // 50 - 70
  VeryLow  = 3   // < 50
};

inline pLDDTCategory categorize_plddt(double score) noexcept {
  if (score > 90.0) return pLDDTCategory::VeryHigh;
  if (score > 70.0) return pLDDTCategory::High;
  if (score > 50.0) return pLDDTCategory::Low;
  return pLDDTCategory::VeryLow;
}

static_assert(static_cast<std::uint8_t>(pLDDTCategory::VeryHigh) == 0, "Enum layout assumed in serialization");
static_assert(sizeof(pLDDTCategory) == sizeof(std::uint8_t), "pLDDTCategory must remain 1 byte");

} // namespace lahuta

#endif // LAHUTA_MODELS_PLDDT_HPP
