#ifndef LAHUTA_GROUPS_HPP
#define LAHUTA_GROUPS_HPP

#include "GraphMol/RWMol.h"
#include "contacts/aromaticity.hpp"
#include "contacts/charges.hpp"
#include "residues.hpp"
#include "entities/records.hpp"
#include <vector>

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
  static std::vector<GroupRec> analyze(const RDKit::RWMol &mol, const Residues &residues);
};

} // namespace lahuta

#endif // LAHUTA_GROUPS_HPP
