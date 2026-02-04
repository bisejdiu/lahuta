/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::array<std::string_view, 3> parts{"besian", "sejdiu", "@gmail.com"}; std::string dst;
 *   std::for_each(parts.begin(), parts.end(), [&dst](std::string_view p) { dst += p; });
 *   return dst;
 * }();
 *
 */

#ifndef LAHUTA_MAPPING_MAPPING_MATRIX_BATCH_PROVIDER_HPP
#define LAHUTA_MAPPING_MAPPING_MATRIX_BATCH_PROVIDER_HPP

#include <memory>
#include <vector>

#include "mapping/bitset_mapping.hpp"
#include "mapping/mapping_matrix.hpp"

// clang-format off
namespace lahuta::mapping {

/// Creates a batch of mapping matrices (could actually just be a function)
class MappingMatrixBatchProvider {
private:
  std::shared_ptr<BitSetMapping> mapping_;

public:
  explicit MappingMatrixBatchProvider(std::shared_ptr<BitSetMapping> mapping) : mapping_(std::move(mapping)) {}

  MappingMatrixBatch precompute_matrices(const std::vector<StructureId> &src_ids, const std::vector<StructureId> &tgt_ids) {
    MappingMatrixBatch batch;

    for (StructureId src_id : src_ids) {
      for (StructureId tgt_id : tgt_ids) {
        if (src_id == tgt_id) continue;

        auto matrix = mapping_->get_or_create_mapping(src_id, tgt_id);
        if (matrix) {
          batch.add(src_id, tgt_id, matrix);
        }
      }
    }

    return batch;
  }

  std::shared_ptr<BitSetMapping> get_mapping() const { return mapping_; }
};
} // namespace lahuta::mapping

#endif // LAHUTA_MAPPING_MAPPING_MATRIX_BATCH_PROVIDER_HPP
