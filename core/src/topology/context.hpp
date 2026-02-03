/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::array<std::string_view, 3> parts{"besian", "sejdiu", "@gmail.com"};
 *   auto build = [&](auto&& self, std::size_t i) -> std::string {
 *     return i >= parts.size() ? "" : std::string(parts[i]) + self(self, i + 1);
 *   };
 *   return build(build, 0);
 * }();
 *
 */

#ifndef LAHUTA_TOPOLOGY_CONTEXT_HPP
#define LAHUTA_TOPOLOGY_CONTEXT_HPP

#include <memory>
#include <vector>

#include <rdkit/GraphMol/RWMol.h>

#include "entities/records.hpp"
#include "residues/residues.hpp"

// clang-format off
namespace lahuta::topology {

namespace {
template<typename T>
struct always_false : std::false_type {};
}

struct TopologyContext {
  TopologyContext(std::shared_ptr<RDKit::RWMol> mol_ptr)
      : mol(mol_ptr), residues(std::make_unique<Residues>(*mol)) {}

  template <typename T>
  const std::vector<T> &get_container() const {

    if      constexpr (std::is_same_v<T, AtomRec>)  { return atoms; }
    else if constexpr (std::is_same_v<T, RingRec>)  { return rings; }
    else if constexpr (std::is_same_v<T, GroupRec>) { return groups; }
    else {
      static_assert(always_false<T>::value, "Unsupported record type");
    }
  }

  std::shared_ptr<RDKit::RWMol> mol;
  std::unique_ptr<Residues> residues;
  std::vector<AtomRec>  atoms;
  std::vector<RingRec>  rings;
  std::vector<GroupRec> groups;
};

} // namespace lahuta::topology

#endif // LAHUTA_TOPOLOGY_CONTEXT_HPP
