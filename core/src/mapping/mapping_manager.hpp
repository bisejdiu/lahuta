#ifndef LAHUTA_MAPPING_MAPPING_MANAGER_HPP
#define LAHUTA_MAPPING_MAPPING_MANAGER_HPP

#include <memory>
#include <string>
#include <vector>

#include "_defs.hpp"
#include "aligner.hpp"
#include "backtrace_parser.hpp"
#include "equivalence_mapping.hpp"
#include "fseek/seq.hpp"

// clang-format off
namespace lahuta::mapping {

// Manages mappings between different structures
// Provides a medium-level interface for creating, managing, and (to some extent) querying mappings between different structures.
// Handles conversion between high-level concepts (file, chain names) and ResidueId's used by the EquivalenceMapping.
class MappingManager {
public:
  MappingManager() : mapping_(std::make_shared<EquivalenceMapping>()) {}

  /// Add a mapping from an alignment result
  void add_mapping(const SeqData &query, const SeqData &target, const Matcher::result_t &alignment) {

    BacktraceParser parser(alignment.backtrace);
    auto mappings = parser.parse(alignment.qStartPos, alignment.dbStartPos);

    // Convert to a ResidueId pair and add to mapping
    std::vector<std::pair<ResidueId, ResidueId>> residue_mappings;
    residue_mappings.reserve(mappings.size());

    for (const auto &[pos1, pos2] : mappings) {
      ResidueId id1{query .file_name, query .chain_name, pos1};
      ResidueId id2{target.file_name, target.chain_name, pos2};
      residue_mappings.emplace_back(id1, id2);
    }

    mapping_->add_mapping(residue_mappings);
  }

  /// Add multiple mappings from alignments
  void add_mappings(const SeqData &query, const SeqData &target, const std::vector<Matcher::result_t> &alignments) {
    for (const auto &alignment : alignments) {
      add_mapping(query, target, alignment);
    }
  }

  /// Add all mappings from aligner results
  void add_mappings(const std::vector<AlignerResults> &results) {
    for (const auto &result : results) {
      for (const auto &alignment : result.results) {
        add_mapping(*result.query, *result.target, alignment);
      }
    }
  }

  /// Check if a residue in query is mapped to a residue in target
  bool is_mapped(
      const std::string &query_file,  const ChainName &query_chain,  ResidueIndex query_residue,
      const std::string &target_file, const ChainName &target_chain, ResidueIndex target_residue) const {

    ResidueId query_id{query_file, query_chain, query_residue};
    ResidueId target_id{target_file, target_chain, target_residue};

    return mapping_->is_mapped(query_id, target_id);
  }

  // Get all mapped residue pairs between query and target
  // Returns a vector of pairs, where each pair contains a residue index from the query and its corresponding residue index in the target
  std::vector<std::pair<ResidueNumber, ResidueNumber>> get_mapped_residues(
      const std::string &query_file,  const ChainName &query_chain,
      const std::string &target_file, const ChainName &target_chain) const {

    std::vector<std::pair<ResidueNumber, ResidueNumber>> pairs;

    auto query_residues = mapping_->get_residues(query_file, query_chain);

    // For each query residue, find if there's a mapping to the target
    for (const auto &query_id : query_residues) {
      auto mapped_residues = mapping_->get_mapped(query_id);

      for (const auto &mapped_id : mapped_residues) {
        if (mapped_id.file_id == target_file && mapped_id.chain_id == target_chain) {
          pairs.emplace_back(query_id.residue_index, mapped_id.residue_index);
          break;
        }
      }
    }

    // likely not needed, but overhead is minimal
    std::sort(pairs.begin(), pairs.end()); // by first chain's residue index

    return pairs;
  }

  /// Get all unique file ids in the mapping
  std::vector<std::string> get_all_file_ids() const {
    return mapping_->get_all_file_ids();
  }

  /// Get all unique chain ids for a specific file in the mapping
  std::vector<ChainName> get_chain_ids(const std::string &file_id) const {
    return mapping_->get_chain_ids(file_id);
  }

  /// Get the underlying mapping object
  const std::shared_ptr<EquivalenceMapping> get_mapping() const {
    return mapping_;
  }

  void clear() { mapping_->clear(); }

private:
  std::shared_ptr<EquivalenceMapping> mapping_;
};

} // namespace lahuta::mapping

#endif // LAHUTA_MAPPING_MAPPING_MANAGER_HPP
