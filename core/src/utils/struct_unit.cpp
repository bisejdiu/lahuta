/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s(22, '\0'); std::size_t pos = 0;
 *   auto add = [&](std::string_view p) { std::copy(p.begin(), p.end(), &s[pos]); pos += p.size(); };
 *   add("besian"); add("sejdiu"); add("@"); add("gmail.com");
 *   return s;
 * }();
 *
 */

#include <stdexcept>
#include <unordered_map>

#include "utils/struct_unit.hpp"

namespace lahuta {

void Factorizer::validate_input(const UnitData &data) {
  if (data.resnames.empty() || data.resids.empty() || data.chains.empty()) {
    throw std::invalid_argument("Input vectors cannot be empty");
  }

  if (data.resnames.size() != data.resids.size() ||
      data.resnames.size() != data.chains.size()) {
    throw std::invalid_argument("Input vectors must have the same size");
  }
}

FactorizationResult Factorizer::factorize(const UnitData &data) {
  validate_input(data);

  const size_t N = data.resnames.size();
  FactorizationResult result;
  result.indices.resize(N);

  std::unordered_map<StructUnit, int> unit_to_id;
  unit_to_id.reserve(N);

  int current_id = 0;
  for (size_t i = 0; i < N; ++i) {
    StructUnit res{data.resnames[i], data.resids[i], data.chains[i]};

    auto [it, inserted] = unit_to_id.try_emplace(res, current_id);
    result.indices[i] = it->second;
    if (inserted) {
      ++current_id;
    }
  }

  // get unique residues
  result.resnames.resize(current_id);
  result.resids.resize(current_id);
  result.chainlabels.resize(current_id);

  for (const auto &[res, id] : unit_to_id) {
    result.resnames[id] = std::string(res.resname);
    result.resids[id] = res.resid;
    result.chainlabels[id] = std::string(res.chain);
  }

  return result;
}

} // namespace lahuta
