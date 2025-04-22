#include "mapping/mapping_matrix.hpp"

namespace lahuta::mapping {

std::vector<ResidueIndex> MappingMatrix::get_source_indices_for_target(ResidueIndex tgt_idx) const {
  std::vector<ResidueIndex> result;
  if (tgt_idx >= tgt_size_) {
    return result;
  }

  // result.reserve(src_size_);
  for (ResidueIndex src_idx = 0; src_idx < src_size_; ++src_idx) {
    if (get(src_idx, tgt_idx)) {
      result.push_back(src_idx);
    }
  }

  return result;
}

std::vector<ResidueIndex> MappingMatrix::get_target_indices_for_source(ResidueIndex src_idx) const {
  std::vector<ResidueIndex> result;
  if (src_idx >= src_size_) {
    return result;
  }

  // result.reserve(tgt_size_);
  for (ResidueIndex tgt_idx = 0; tgt_idx < tgt_size_; ++tgt_idx) {
    if (get(src_idx, tgt_idx)) {
      result.push_back(tgt_idx);
    }
  }

  return result;
}

std::vector<std::pair<ResidueIndex, ResidueIndex>> MappingMatrix::get_all_mappings() const {
  std::vector<std::pair<ResidueIndex, ResidueIndex>> result;
  result.reserve(matrix_.count());

  for (ResidueIndex src_idx = 0; src_idx < src_size_; ++src_idx) {
    for (ResidueIndex tgt_idx = 0; tgt_idx < tgt_size_; ++tgt_idx) {
      if (get(src_idx, tgt_idx)) {
        result.emplace_back(src_idx, tgt_idx);
      }
    }
  }

  return result;
}

} // namespace lahuta::mapping
