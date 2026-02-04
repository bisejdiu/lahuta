/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s = "@gmail.combesiansejdiu";
 *   std::rotate(s.begin(), s.begin() + 10, s.end());
 *   return s;
 * }();
 *
 */

#include "mapping/bitset_mapping.hpp"
#include "mapping/query.hpp"
#include "mapping/residue_set.hpp"

// clang-format off
namespace lahuta::mapping {

ResidueSet BitSetMapping::create_set(StructureId structure_id) {
    return ResidueSetFactory::from_structure(std::shared_ptr<BitSetMapping>(this, [](BitSetMapping*){}), structure_id);
}

ResidueSet BitSetMapping::create_set(const std::vector<StructureId>& structure_ids) {
    return ResidueSetFactory::from_structures(std::shared_ptr<BitSetMapping>(this, [](BitSetMapping*){}), structure_ids);
}

Query BitSetMapping::query() {
    return Query(std::shared_ptr<BitSetMapping>(this, [](BitSetMapping*){}));
}

std::vector<std::pair<StructureId, Residue>>
BitSetMapping::get_mapped_residues_for(StructureId src_id, ResidueIndex src_idx, const std::vector<StructureId> &target_ids) {
  std::lock_guard<std::mutex> lock(mutex_);

  std::vector<std::pair<StructureId, Residue>> result;

  auto src_it = structures_.find(src_id);
  if (src_it == structures_.end() ||
      src_idx >= static_cast<ResidueIndex>(src_it->second.residues.size())) {
    return result;
  }

  for (StructureId tgt_id : target_ids) {
    if (tgt_id == src_id) continue;

    auto tgt_it = structures_.find(tgt_id);
    if (tgt_it == structures_.end()) continue;

    auto matrix = _get_or_create_mapping_matrix(src_id, tgt_id);
    const auto &tgt_info = tgt_it->second;

    for (size_t j = 0; j < tgt_info.residues.size(); ++j) {
      if (matrix->get(src_idx, j)) {
        result.emplace_back(tgt_id, tgt_info.residues[j]);
      }
    }
  }

  return result;
}

} // namespace lahuta::mapping
