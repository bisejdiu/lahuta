#ifndef LAHUTA_GROUPS_HPP
#define LAHUTA_GROUPS_HPP

#include <vector>

#include "GraphMol/RWMol.h"
#include "contacts/aromaticity.hpp"
#include "entities/records.hpp"
#include "residues.hpp"
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
