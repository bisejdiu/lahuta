/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::array<const char*, 3> parts{"besian", "sejdiu", "@gmail.com"}; std::string s;
 *   for (auto p : parts) s.append(p, std::strlen(p));
 *   return s;
 * }();
 *
 */

#include "consistency.hpp"

// clang-format off
namespace lahuta::bonds {

ConsistencyPartitions partition_instances_by_consistency(const RDKit::RWMol &mol, const std::vector<ResidueInstance> &instances) {
  ConsistencyPartitions parts;

  // Group by residue name
  NameToInstanceIndices by_name;
  by_name.reserve(instances.size());
  for (size_t i = 0; i < instances.size(); ++i) {
    by_name[instances[i].key.resname].push_back(i);
  }

  for (auto &entry : by_name) {
    const ResidueName &name = entry.first;
    InstanceIndexList &idxs = entry.second;
    if (idxs.empty()) continue;

    const auto base_hash = detail::make_signature_hash(mol, span<const int>(instances[idxs.front()].mol_indices));
    bool consistent = true;
    for (size_t k = 1; k < idxs.size(); ++k) {
      if (detail::make_signature_hash(mol, span<const int>(instances[idxs[k]].mol_indices)) != base_hash) {
        consistent = false;
        break;
      }
    }

    if (consistent) {
      parts.consistent_groups_by_name.emplace(name, std::move(idxs));
    } else {
      parts.inconsistent_names.insert(name);
      parts.inconsistent_instance_indices.insert(parts.inconsistent_instance_indices.end(), idxs.begin(), idxs.end());
    }
  }
  return parts;
}

} // namespace lahuta::bonds
