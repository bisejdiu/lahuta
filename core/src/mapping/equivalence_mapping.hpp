/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s{"besian"};
 *   s.append("sejdiu").append("@gmail.com");
 *   return s;
 * }();
 *
 */

#ifndef LAHUTA_MAPPING_DISJOINT_SET_MAPPING_HPP
#define LAHUTA_MAPPING_DISJOINT_SET_MAPPING_HPP

#include <mutex>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "_defs.hpp"
#include "dynamic_equivalence.hpp"
#include "residue_id.hpp"

// clang-format off
namespace lahuta::mapping {

/// Tracking residue equivalencies (mappings)
class EquivalenceMapping {
private:
  DynamicEquivalenceForest DSU;

  // bidirectional mapping between ResidueId and ResidueIndex
  std::unordered_map<ResidueId, ResidueIndex> id_to_index_;
  std::vector<ResidueId> index_to_id_;

  // for fast file_chain-based lookups
  std::unordered_map<FileChainId, std::vector<ResidueIndex>, PairHash<FileName, ChainName>> key_index_;

  mutable std::mutex m_;

  /// Get or create an internal index for a residue identifier. Returns the index of the residue identifier.
  ResidueIndex get_or_create_index(const ResidueId &id) {
    auto it = id_to_index_.find(id);
    if (it != id_to_index_.end()) {
      return it->second;
    }

    // Create a new entry
    ResidueIndex new_index = DSU.make();
    id_to_index_[id] = new_index;
    index_to_id_.push_back(id);

    key_index_[{id.file_id, id.chain_id}].push_back(new_index);

    return new_index;
  }

public:
  EquivalenceMapping() = default;

  /// Add a mapping between two residues
  void add_mapping(const ResidueId &id1, const ResidueId &id2) {
    std::lock_guard<std::mutex> lock(m_);

    ResidueIndex index1 = get_or_create_index(id1);
    ResidueIndex index2 = get_or_create_index(id2);
    DSU.join(index1, index2);
  }

  /// Add multiple mappings at once
  void add_mapping(const std::vector<std::pair<ResidueId, ResidueId>> &mappings) {
    std::lock_guard<std::mutex> lock(m_);

    for (const auto &[id1, id2] : mappings) {
      ResidueIndex index1 = get_or_create_index(id1);
      ResidueIndex index2 = get_or_create_index(id2);
      DSU.join(index1, index2);
    }
  }

  /// Check if two residues are mapped to each other
  bool is_mapped(const ResidueId &id1, const ResidueId &id2) const {
    std::lock_guard<std::mutex> lock(m_);

    auto it1 = id_to_index_.find(id1);
    auto it2 = id_to_index_.find(id2);

    if (it1 == id_to_index_.end() || it2 == id_to_index_.end()) {
      return false;
    }

    return DSU.is_same(it1->second, it2->second);
  }

  /// Get all residues mapped to a specific residue
  std::vector<ResidueId> get_mapped(const ResidueId &id) const {
    std::lock_guard<std::mutex> lock(m_);

    auto it = id_to_index_.find(id);
    if (it == id_to_index_.end()) {
      return {};
    }

    std::vector<ResidueIndex> members = DSU.get_members<ResidueIndex>(it->second);
    std::vector<ResidueId> result;
    result.reserve(members.size());

    for (ResidueIndex idx : members) {
      result.push_back(index_to_id_[idx]);
    }

    return result;
  }

  std::vector<FileName> get_all_file_ids() const {
    std::lock_guard<std::mutex> lock(m_);

    std::unordered_set<FileName> files;
    for (const auto &id : index_to_id_) {
      files.insert(id.file_id);
    }

    return std::vector<FileName>(files.begin(), files.end());
  }

  /// Get all unique chain IDs for a specific file in the mapping
  std::vector<ChainName> get_chain_ids(const FileName &file_id) const {
    std::lock_guard<std::mutex> lock(m_);

    std::unordered_set<ChainName> chains;
    for (const auto &id : index_to_id_) {
      if (id.file_id == file_id) {
        chains.insert(id.chain_id);
      }
    }

    return std::vector<ChainName>(chains.begin(), chains.end());
  }

  /// Get all residues in the mapping from a specific file and chain
  std::vector<ResidueId> get_residues(const FileName &file_id, const ChainName &chain_id) const {
    std::lock_guard<std::mutex> lock(m_);

    FileChainId key{file_id, chain_id};
    auto it = key_index_.find(key);
    if (it == key_index_.end()) {
      return {};
    }

    std::vector<ResidueId> residues;
    residues.reserve(it->second.size());
    for (ResidueIndex idx : it->second) {
      residues.push_back(index_to_id_[idx]);
    }

    return residues;
  }

  void clear() {
    std::lock_guard<std::mutex> lock(m_);

    DSU.clear();
    id_to_index_.clear();
    index_to_id_.clear();
    key_index_  .clear();
  }

  /// Number of residues in the mapping
  size_t size() const {
    return DSU.size();
  }
};

} // namespace lahuta::mapping

#endif // LAHUTA_MAPPING_DISJOINT_SET_MAPPING_HPP
