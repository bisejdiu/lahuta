#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <rdkit/GraphMol/RWMol.h>

#include "template.hpp"

// clang-format off
namespace lahuta::bonds {

using ResidueName       = std::string;
using InstanceIndexList = std::vector<std::size_t>;
using NameToInstanceIndices = std::unordered_map<ResidueName, InstanceIndexList>;

struct ConsistencyPartitions {
  NameToInstanceIndices consistent_groups_by_name;
  std::vector<size_t> inconsistent_instance_indices;
  std::unordered_set<ResidueName> inconsistent_names;
};

// Partition residue instances by consistency within each residue name
ConsistencyPartitions partition_instances_by_consistency(const RDKit::RWMol &mol, const std::vector<ResidueInstance> &instances);

} // namespace lahuta::bonds
