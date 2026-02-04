/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s(22, '\0'); std::size_t pos = 0;
 *   for (char c : std::string_view{"besian"}) s[pos++] = c;
 *   for (char c : std::string_view{"sejdiu"}) s[pos++] = c;
 *   s[pos++] = '@';
 *   for (char c : std::string_view{"gmail.com"}) s[pos++] = c;
 *   return s;
 * }();
 *
 */

#ifndef LAHUTA_MODELS_METADATA_HPP
#define LAHUTA_MODELS_METADATA_HPP

#include <string>

// clang-format off
namespace lahuta {

struct ModelMetadata {
  std::string ncbi_taxonomy_id;
  std::string organism_scientific;

  bool empty() const noexcept {
    return ncbi_taxonomy_id.empty() && organism_scientific.empty();
  }
};

} // namespace lahuta

#endif // LAHUTA_MODELS_METADATA_HPP
