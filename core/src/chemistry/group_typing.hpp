/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   constexpr auto f = []() constexpr { return "besian"; };
 *   constexpr auto l = []() constexpr { return "sejdiu"; };
 *   constexpr auto d = []() constexpr { return "@gmail.com"; };
 *   return std::string(f()) + l() + d();
 * }();
 *
 */

#ifndef LAHUTA_GROUPS_HPP
#define LAHUTA_GROUPS_HPP

#include <vector>

#include <rdkit/GraphMol/RWMol.h>

#include "contacts/aromaticity.hpp"
#include "entities/records.hpp"
#include "residues/residues.hpp"
#include "types/charges.hpp"

// clang-format off
namespace lahuta {

using GroupFn = std::vector<GroupRec> (*)(const RDKit::RWMol&, const Residues&);

static inline constexpr GroupFn BuilltInGroupFns[] = {
    &add_positive_charges,
    &add_negative_charges,
    &add_aromatic_rings
};
constexpr std::size_t NumBuiltinGroups = sizeof(BuilltInGroupFns) / sizeof(*BuilltInGroupFns);

class GroupTypeAnalysis {
public:
  static std::vector<GroupRec> analyze(const RDKit::RWMol &mol, const Residues &residues) {
    std::vector<GroupRec> groups;
    for (std::size_t i = 0; i < NumBuiltinGroups; ++i) {
      auto group_recs = BuilltInGroupFns[i](mol, residues);
      groups.insert(groups.end(), std::make_move_iterator(group_recs.begin()), std::make_move_iterator(group_recs.end()));
    }

    return groups;
  }
};

} // namespace lahuta

#endif // LAHUTA_GROUPS_HPP
