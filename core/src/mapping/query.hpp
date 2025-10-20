#ifndef LAHUTA_MAPPING_QUERY_HPP
#define LAHUTA_MAPPING_QUERY_HPP

#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "_defs.hpp"
#include "bitset_mapping.hpp"
#include "residue_set.hpp"

// clang-format off
namespace lahuta::mapping {

/// Query interface for residue mappings
class Query {
public:
  explicit Query(std::shared_ptr<BitSetMapping> mapping) : mapping_(std::move(mapping)) {}

  explicit Query(ResidueSet residue_set)
      : mapping_(residue_set.getMapping()), initial_set_(std::move(residue_set)) {}

  /// Set the source structure(s) to query from
  Query &from(StructureId structure_id) {
    from_structures_.clear();
    from_structures_.push_back(structure_id);
    return *this;
  }

  Query &from(const std::vector<StructureId> &structure_ids) {
    from_structures_ = structure_ids;
    return *this;
  }

  Query &from(const ResidueSet &residue_set) {
    // Ensure both are using the same mapping
    if (mapping_ != residue_set.getMapping()) {
      throw std::invalid_argument("Cannot query from a ResidueSet with a different mapping");
    }

    initial_set_ = residue_set;
    return *this;
  }

  /// Set the target structure(s) to find mappings to
  Query &to(StructureId structure_id) {
    to_structures_.clear();
    to_structures_.push_back(structure_id);
    return *this;
  }

  Query &to(const std::vector<StructureId> &structure_ids) {
    to_structures_ = structure_ids;
    return *this;
  }

  // Filter residues based on a predicate
  template <typename Predicate>
  Query &where(Predicate &&predicate) {
    static_assert(std::is_invocable_r_v<bool, Predicate, StructureId, ResidueIndex>,
                  "Predicate must be callable as bool(StructureId,ResidueIndex)");
    predicates_.push_back(std::forward<Predicate>(predicate));
    return *this;
  }

  /// Filter residues based on the residue name
  Query &where_residue_name(const ResidueName &name) {
    return where([this, name](StructureId sid, ResidueIndex idx) {
      auto res = mapping_->get_residue(sid, idx);
      return res && res->name == name;
    });
  }

  /// Filter residues based on a range of residue numbers
  Query &where_residue_number_in_range(ResidueNumber start, ResidueNumber end) {
    return where([this, start, end](StructureId sid, ResidueIndex idx) {
      auto res = mapping_->get_residue(sid, idx);
      return res && res->number >= start && res->number <= end;
    });
  }

  /// Filter residues based on a specific chain
  Query &where_chain(const std::string &chain) {
    return where([this, chain](StructureId sid, ResidueIndex idx) {
      auto res = mapping_->get_residue(sid, idx);
      return res && res->chain_id == chain;
    });
  }

  /// Execute the query and return a set of residues matching the specified criteria
  ///
  /// This is the method that applies all query filters (from/to/where clauses)
  /// and returns a ResidueSet that can be used for further set operations.
  /// All select methods internally call this method first.
  ResidueSet execute() {
    // use the initial set if we have one, otherwise create a new one
    ResidueSet result_set = initial_set_.has_value()
                                ? initial_set_.value()
                                : ResidueSetFactory::from_structures(mapping_, from_structures_);

    // If targeting specific structures, get mapped residues
    if (!to_structures_.empty()) {
      result_set = get_mapped_residues(result_set, to_structures_);
    }

    // Apply all predicates
    for (const auto &predicate : predicates_) {
      result_set = result_set.filter(predicate);
    }

    return result_set;
  }

  /// Transform each matched residue using a custom transformation function
  template <typename T, typename TransformFunc>
  std::vector<T> select(TransformFunc &&transform) {
    static_assert(std::is_invocable_r_v<T, TransformFunc, StructureId, ResidueIndex>,
                  "TransformFunc must be callable as T(StructureId,ResidueIndex)");

    auto result_set = execute();
    return result_set.map<T>(std::forward<TransformFunc>(transform)); // Map each residue to the desired property
  }

  /// A specialization of select<T>() returning the structure ID and residue index for each matched residue.
  std::vector<std::pair<StructureId, ResidueIndex>> select_indices() {
    return select<std::pair<StructureId, ResidueIndex>>(
        [](StructureId sid, ResidueIndex idx) { return std::make_pair(sid, idx); });
  }

  /// A specialization of select<T>() returning the structure ID and Residue object for each matched residue.
  std::vector<std::pair<StructureId, Residue>> select_residues() {
    return select<std::pair<StructureId, Residue>>(
        [this](StructureId sid, ResidueIndex idx) {
          auto res = mapping_->get_residue(sid, idx);
          return std::make_pair(sid, res.value_or(Residue{}));
        });
  }

private:
  /// Get residues from source_set that map to any structure in target_ids
  ResidueSet get_mapped_residues(const ResidueSet &source_set, const std::vector<StructureId> &target_ids) {
    // For each structure in the source set, find residues that map to any target
    std::vector<ResidueSet::StructureResidues> mapped_collections;

    for (const auto &collection : source_set.get_collection()) {
      StructureId src_id = collection.structure_id;
      std::vector<ResidueIndex> mapped_indices;

      for (ResidueIndex src_idx : collection.indices) {

        // find if this residue is mapped to any of the target structures
        bool is_mapped = false;
        for (StructureId tgt_id : target_ids) {
          if (mapping_->is_residue_mapped(src_id, src_idx, tgt_id)) {
            is_mapped = true;
            break;
          }
        }

        if (is_mapped) {
          mapped_indices.push_back(src_idx);
        }
      }

      if (!mapped_indices.empty()) {
        mapped_collections.emplace_back(src_id, std::move(mapped_indices));
      }
    }

    return ResidueSet(mapping_, std::move(mapped_collections));
  }

  std::optional<ResidueSet>       initial_set_;       // Optional initial set of residues
  std::vector<StructureId>        from_structures_;   // Source structures to query from
  std::vector<StructureId>        to_structures_;     // Target structures to find mappings to
  std::shared_ptr<BitSetMapping>  mapping_;
  std::vector<std::function<bool(StructureId, ResidueIndex)>> predicates_; // List of predicates to filter residues
};

} // namespace lahuta::mapping

#endif // LAHUTA_MAPPING_QUERY_HPP 
