#ifndef LAHUTA_MAPPING_MAPPING_MATRIX_HPP
#define LAHUTA_MAPPING_MAPPING_MATRIX_HPP

#include "_defs.hpp"
#include <boost/dynamic_bitset.hpp>
#include <cstdint>
#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>

// clang-format off
namespace lahuta::mapping {

/// Handles a mapping between residues of two structures
class MappingMatrix {
private:
  StructureId src_id_;
  StructureId tgt_id_;
  ResidueIndex src_size_;
  ResidueIndex tgt_size_;
  boost::dynamic_bitset<> matrix_;

public:
  MappingMatrix(StructureId src, StructureId tgt, ResidueIndex src_size, ResidueIndex tgt_size)
      : src_id_(src), tgt_id_(tgt), matrix_(src_size * tgt_size), src_size_(src_size), tgt_size_(tgt_size) {}

  const StructureId  get_src_id()   const noexcept { return src_id_; }
  const StructureId  get_tgt_id()   const noexcept { return tgt_id_; }
  const ResidueIndex get_src_size() const noexcept { return src_size_; }
  const ResidueIndex get_tgt_size() const noexcept { return tgt_size_; }

  /// Get the mapping value at specified position
  inline bool get(ResidueIndex src_idx, ResidueIndex tgt_idx) const noexcept {
    ResidueIndex flat_idx = src_idx * tgt_size_ + tgt_idx;
    return matrix_[flat_idx];
  }

  /// Set the mapping value at specified position
  void set(ResidueIndex src_idx, ResidueIndex tgt_idx, bool value = true) {
    ResidueIndex flat_idx = src_idx * tgt_size_ + tgt_idx;
    matrix_[flat_idx] = value;
  }

  /// Number of mapped pairs in the matrix
  auto count() const { return matrix_.count(); }

  /// Does the matrix contain any mappings?
  bool has_mappings() const { return matrix_.any(); }

  /// Get all source indices that map to the given target index
  std::vector<ResidueIndex> get_source_indices_for_target(ResidueIndex tgt_idx) const;

  /// Get all target indices that map to the given source index
  std::vector<ResidueIndex> get_target_indices_for_source(ResidueIndex src_idx) const;

  /// Get all pairs of indices that are mapped
  std::vector<std::pair<ResidueIndex, ResidueIndex>> get_all_mappings() const;

  /// Create a unique key for storing mapping matrices
  static uint64_t create_matrix_key(StructureId src_id, StructureId tgt_id) {
    return (static_cast<uint64_t>(src_id) << 32) | static_cast<uint64_t>(tgt_id);
  }
};

// TODO: Further (maybe micro) optimizations:
//  1. Remove pair construction using heterogeneous lookup (does not seem to be much faster)
//  2. Remove atomic ref counting for shared_ptrs (use raw or unique_ptrs) (don't see it making a difference)
//  3. Use alternative implementations of std::unordered_map (likely to provide noticeable performance improvements)

/// Handles a batch of mapping matrices between multiple structures
class MappingMatrixBatch {
private:
  using StructureIdPair = std::pair<StructureId, StructureId>;
  using SPairHash       = PairHash <StructureId, StructureId>;

  std::unordered_map<StructureIdPair, std::shared_ptr<MappingMatrix>, SPairHash> matrices_;

public:
  MappingMatrixBatch() = default;

  /// if we already have matrices
  explicit MappingMatrixBatch(std::unordered_map<StructureIdPair, std::shared_ptr<MappingMatrix>, SPairHash> matrices)
      : matrices_(std::move(matrices)) {}

  /// Get a mapping matrix for a specific pair of structures
  [[gnu::always_inline, gnu::hot, nodiscard]]
  inline auto get(StructureId src_id, StructureId tgt_id) const noexcept {
    auto it = matrices_.find(std::make_pair(src_id, tgt_id));
    if (it != matrices_.end()) {
      return it->second;
    }
    return std::shared_ptr<MappingMatrix>(nullptr);
  }

  /// Do we have a mapping matrix for the given pair of structures?
  bool contains(StructureId src_id, StructureId tgt_id) const {
    return matrices_.find(std::make_pair(src_id, tgt_id)) != matrices_.end();
  }

  /// Add a mapping matrix to the batch
  void add(StructureId src_id, StructureId tgt_id, std::shared_ptr<MappingMatrix> matrix) {
    matrices_[std::make_pair(src_id, tgt_id)] = std::move(matrix);
  }

  // Add all matrices from another batch
  void merge(const MappingMatrixBatch &other) {
    for (const auto &[pair, matrix] : other.matrices_) {
      matrices_[pair] = matrix;
    }
  }
};

} // namespace lahuta::mapping

#endif // LAHUTA_MAPPING_MAPPING_MATRIX_HPP
