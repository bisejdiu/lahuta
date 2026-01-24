#ifndef LAHUTA_MAPPING_BITSET_MAPPING_HPP
#define LAHUTA_MAPPING_BITSET_MAPPING_HPP

#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include <boost/dynamic_bitset.hpp>

#include "_defs.hpp"
#include "mapping_manager.hpp"
#include "mapping_matrix.hpp"
#include "residues/residues.hpp"
#include "topology.hpp"

// NOTE: We make copies of the residues in the mapping. Residues are POD objects,
// so it's not clear to me if the tradoff of the better performance and memory footprint,
// compared to the hit on code readability is worth it.

// clang-format off
namespace lahuta::mapping {
class Query;
class ResidueSet;

/// Efficient residue mapping using bit vectors for set operations
class BitSetMapping {
public:
  struct StructureInfo { // A single structure with its residues
    StructureId id;
    std::string file_name;
    std::string chain_name;
    const Topology *topology;
    std::vector<Residue> residues;

    // NOTE: memory overhead and actual real-world benefits have not been measured
    std::map<ChainResidue, ResidueIndex, ChainResidueCompare> chain_residue_index_map; // for faster residue access
  };

private:
  mutable std::mutex mutex_;
  StructureId next_id_ = 1;                 // global residue counter for unique ids
  std::shared_ptr<MappingManager> manager_; // mapping impl

  std::unordered_map<StructureId, StructureInfo> structures_;                     // structure registry
  std::unordered_map<uint64_t, std::shared_ptr<MappingMatrix>> mapping_matrices_; // matrix mapping cache

  /// Get or create a mapping matrix between two structures. Caller should handle locking.
  std::shared_ptr<MappingMatrix> _get_or_create_mapping_matrix(StructureId src_id, StructureId tgt_id) {
    uint64_t key = MappingMatrix::create_matrix_key(src_id, tgt_id);
    auto it = mapping_matrices_.find(key);
    if (it != mapping_matrices_.end()) {
      return it->second;
    }

    const auto &src_info = structures_[src_id];
    const auto &tgt_info = structures_[tgt_id];

    auto matrix = std::make_shared<MappingMatrix>( // create new mapping matrix
        src_id,
        tgt_id,
        src_info.residues.size(),
        tgt_info.residues.size());

    auto pairs = manager_->get_mapped_residues( // get the mapped residues to populate the matrix
        src_info.file_name,
        src_info.chain_name,
        tgt_info.file_name,
        tgt_info.chain_name);

    for (const auto &[src_idx, tgt_idx] : pairs) {
      matrix->set(src_idx, tgt_idx);
    }

    mapping_matrices_[key] = matrix;
    return matrix;
  }

public:
  BitSetMapping() : manager_(std::make_shared<MappingManager>()) {}

  explicit BitSetMapping(std::shared_ptr<MappingManager> manager) : manager_(std::move(manager)) {}

  /// Register a structure; returns its unique ID
  StructureId register_structure(const std::string &file_name, const std::string &chain_name, const Topology &topology) {
    std::lock_guard<std::mutex> lock(mutex_);

    StructureId id = next_id_++;
    StructureInfo info;
    info.file_name  = file_name;
    info.chain_name = chain_name;
    info.topology   = &topology;
    info.id         = id;

    // Get residues and create global-to-local mapping
    const auto &residues = topology.get_residues().get_residues();
    info.residues.reserve(residues.size());

    for (u32 i = 0; i < static_cast<u32>(residues.size()); ++i) {
      info.residues.push_back(residues[i]);

      // Populate lookup maps
      const auto& residue = residues[i];
      auto key = std::make_tuple(residue.chain_id, residue.number, residue.name);
      info.chain_residue_index_map[key] = i;
    }

    structures_[info.id] = std::move(info);
    return id;
  }

  /// Get a mapping matrix between two structures
  std::shared_ptr<MappingMatrix> get_or_create_mapping(StructureId src_id, StructureId tgt_id) {
    std::lock_guard<std::mutex> lock(mutex_);
    return _get_or_create_mapping_matrix(src_id, tgt_id);
  }

  /// Get all registered structure ids
  std::vector<StructureId> get_all_structure_ids() const {
    std::lock_guard<std::mutex> lock(mutex_);
    std::vector<StructureId> ids;
    ids.reserve(structures_.size());

    for (const auto& [id, _] : structures_) {
      ids.push_back(id);
    }

    return ids;
  }

  std::optional<StructureInfo> get_structure_info(StructureId id) const {
    std::lock_guard<std::mutex> lock(mutex_);
    auto it = structures_.find(id);
    if (it != structures_.end()) {
      return it->second;
    }
    return std::nullopt;
  }

  /// Get the indices of all residues in a structure
  std::vector<ResidueIndex> get_residue_indices(StructureId id) const {
    std::lock_guard<std::mutex> lock(mutex_);

    auto it = structures_.find(id);
    if (it == structures_.end()) {
      return {};
    }

    const auto& info = it->second;
    std::vector<ResidueIndex> indices(info.residues.size());
    for (size_t i = 0; i < info.residues.size(); ++i) {
      indices[i] = static_cast<ResidueIndex>(i);
    }

    return indices;
  }

  /// Get a specific residue by its mapping index within a structure
  std::optional<Residue> get_residue(StructureId id, ResidueIndex index) const {
    std::lock_guard<std::mutex> lock(mutex_);
    auto it = structures_.find(id);
    if (it == structures_.end() || index >= static_cast<ResidueIndex>(it->second.residues.size())) {
      return std::nullopt;
    }

    return it->second.residues[index];
  }

  /// Check if a residue in the source structure is mapped to any residue in the target structure
  bool is_residue_mapped(StructureId src_id, ResidueIndex src_idx, StructureId tgt_id) {
    std::lock_guard<std::mutex> lock(mutex_);

    auto src_it = structures_.find(src_id);
    auto tgt_it = structures_.find(tgt_id);

    if (
      src_it == structures_.end() ||
      tgt_it == structures_.end() ||
      src_idx >= static_cast<ResidueIndex>(src_it->second.residues.size())) {
      return false;
    }

    uint64_t key = MappingMatrix::create_matrix_key(src_id, tgt_id);
    auto matrix_it = mapping_matrices_.find(key);

    if (matrix_it == mapping_matrices_.end()) {
      auto matrix = _get_or_create_mapping_matrix(src_id, tgt_id);

      // Check if src_idx is mapped to any residue in tgt_id
      for (size_t j = 0; j < tgt_it->second.residues.size(); ++j) {
        if (matrix->get(src_idx, j)) {
          return true;
        }
      }
    } else {
      // we have the matrix
      auto matrix = matrix_it->second;

      // Check if src_idx is mapped to any residue in tgt_id
      for (size_t j = 0; j < tgt_it->second.residues.size(); ++j) {
        if (matrix->get(src_idx, j)) {
          return true;
        }
      }
    }

    return false;
  }

  /// Check if a residue in the source structure is mapped to a residue in the target structure
  bool are_residues_mapped(StructureId src_id, ResidueIndex src_idx, StructureId tgt_id, ResidueIndex tgt_idx) {
    std::lock_guard<std::mutex> lock(mutex_);

    auto src_it = structures_.find(src_id);
    auto tgt_it = structures_.find(tgt_id);

    if (src_it == structures_.end() || tgt_it == structures_.end() || 
        src_idx >= static_cast<ResidueIndex>(src_it->second.residues.size()) ||
        tgt_idx >= static_cast<ResidueIndex>(tgt_it->second.residues.size())) {
      return false;
    }

    auto matrix = _get_or_create_mapping_matrix(src_id, tgt_id);
    return matrix->get(src_idx, tgt_idx);
  }

  /// Check if a residue in the source structure is mapped to a residue in the target structure
  bool are_residues_mapped(StructureId src_id, const Residue &src_residue, StructureId tgt_id, const Residue &tgt_residue) {
    std::lock_guard<std::mutex> lock(mutex_);

    ResidueIndex src_idx = find_residue_index(src_id, src_residue.chain_id, src_residue.number, src_residue.name);
    if (src_idx == static_cast<ResidueIndex>(-1)) return false;

    ResidueIndex tgt_idx = find_residue_index(tgt_id, tgt_residue.chain_id, tgt_residue.number, tgt_residue.name);
    if (tgt_idx == static_cast<ResidueIndex>(-1)) return false;

    // Check the mapping matrix
    auto matrix = _get_or_create_mapping_matrix(src_id, tgt_id);
    return matrix->get(src_idx, tgt_idx);
  }


  /// Get all residues that are mapped to a specific residue in the target structure
  std::vector<std::pair<Residue, Residue>> get_mapped_residues(StructureId src_id, StructureId tgt_id) {
    std::lock_guard<std::mutex> lock(mutex_);

    auto matrix = _get_or_create_mapping_matrix(src_id, tgt_id);
    const auto &src_info = structures_.at(src_id);
    const auto &tgt_info = structures_.at(tgt_id);

    std::vector<std::pair<Residue, Residue>> result;

    // find mapped pairs
    for (size_t i = 0; i < src_info.residues.size(); ++i) {
      for (size_t j = 0; j < tgt_info.residues.size(); ++j) {
        if (matrix->get(i, j)) {
          result.emplace_back(src_info.residues[i], tgt_info.residues[j]);
        }
      }
    }

    return result;
  }

  // Get all residues in target_ids that map to the specified source residue
  std::vector<std::pair<StructureId, Residue>>
  get_mapped_residues_for(StructureId src_id, ResidueIndex src_idx, const std::vector<StructureId>& target_ids);

  /// Fast lookup of residue index by chain_id, number, and name
  ResidueIndex find_residue_index(StructureId id, const std::string& chain_id, ResidueNumber number, const ResidueName& name) const {
    std::lock_guard<std::mutex> lock(mutex_);
    auto it = structures_.find(id);
    if (it == structures_.end()) {
      return static_cast<ResidueIndex>(-1);
    }

    const auto& info = it->second;
    auto key = std::make_tuple(chain_id, number, name);
    auto map_it = info.chain_residue_index_map.find(key);
    if (map_it != info.chain_residue_index_map.end()) {
      return map_it->second;
    }

    return static_cast<ResidueIndex>(-1);
  }

  /// Create a ResidueSet from a single structure
  ResidueSet create_set(StructureId structure_id);

  /// Create a ResidueSet from multiple structures
  ResidueSet create_set(const std::vector<StructureId>& structure_ids);

  /// Create a Query object
  Query query();

  /// Get the underlying MappingManager
  std::shared_ptr<MappingManager> get_manager() const { return manager_; }

  /// Add a single alignment result to the mapping
  void add_alignment_result(const AlignerResults& result) {
    std::lock_guard<std::mutex> lock(mutex_);
    // Process all alignments in the result, not just the first one
    for (const auto& alignment : result.results) {
      manager_->add_mapping(*result.query, *result.target, alignment);
    }
  }

  void add_alignment_result(const AlignerResultsX& result) {
    std::lock_guard<std::mutex> lock(mutex_);
    // Process all alignments in the result, not just the first one
    for (const auto& alignment : result.results) {
      manager_->add_mapping(*result.query, *result.target, alignment);
    }
  }
};

// Create a BitSetMapping from alignment results
inline std::shared_ptr<BitSetMapping> create_bit_set_mapping(const std::vector<AlignerResults> &results) {
  auto manager = std::make_shared<MappingManager>();
  manager->add_mappings(results);
  return std::make_shared<BitSetMapping>(manager);
}

} // namespace lahuta::mapping

#endif // LAHUTA_MAPPING_BITSET_MAPPING_HPP
