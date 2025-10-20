#include "mapping/bitset_mapping.hpp"
#include "mapping/mapping_matrix.hpp"
#include "mapping/mapping_matrix_batch_provider.hpp"
#include "mapping/residue_set.hpp"

// clang-format off
namespace lahuta::mapping {

namespace detail {

/// Extract structure IDs from a collection of StructureResidues
inline std::vector<StructureId> extract_structure_ids(const std::vector<ResidueSet::StructureResidues> &collections) {
  std::vector<StructureId> result;
  result.reserve(collections.size());
  for (const auto &collection : collections) {
    result.push_back(collection.structure_id);
  }

  return result;
}

/// Check if a residue exists in a collection for the same structure
inline bool exists_in_same_structure(const ResidueSet::StructureResidues &other_collection, StructureId our_sid, ResidueIndex our_idx) {
  if (other_collection.structure_id != our_sid) {
    return false;
  }

  return std::binary_search(other_collection.indices.begin(), other_collection.indices.end(), our_idx);
}

/// Check if a residue maps to any residue in another collection using direct mapping
inline bool maps_to_any_in_collection(
    const std::shared_ptr<BitSetMapping> &mapping, StructureId our_sid, ResidueIndex our_idx,
    const ResidueSet::StructureResidues &other_collection) {
  StructureId other_sid = other_collection.structure_id;

  if (our_sid == other_sid) return false; // NOTE: we handle this case separately

  for (ResidueIndex other_idx : other_collection.indices) {
    if (mapping->are_residues_mapped(our_sid, our_idx, other_sid, other_idx)) {
      return true;
    }
  }

  return false;
}

/// Check if a residue maps to any residue in another collection using pre-computed matrix
inline bool maps_to_any_in_collection(
    const std::shared_ptr<MappingMatrix> &matrix, ResidueIndex our_idx,
    const ResidueSet::StructureResidues &other_collection) {
  if (!matrix) return false;

  for (ResidueIndex other_idx : other_collection.indices) {
    if (matrix->get(our_idx, other_idx)) {
      return true;
    }
  }

  return false;
}

/// Verify that both sets use the same underlying mapping
inline void verify_mappings(
    const std::shared_ptr<BitSetMapping> &first, const std::shared_ptr<BitSetMapping> &second,
    const char *operation_name) {
  if (first != second) {
    throw std::invalid_argument(
        std::string("Cannot perform ") + operation_name + " on ResidueSet objects with different mappings");
  }
}

/// abstracts the common logic between intersection and difference
template <typename InclusionPredicate>
inline ResidueSet rs_set_operation(const ResidueSet &first, const ResidueSet &second, InclusionPredicate &&predicate) {
  verify_mappings(first.getMapping(), second.getMapping(), "set operation");

  std::vector<ResidueSet::StructureResidues> result_collections;

  // gor each structure in our set
  for (const auto &our_collection : first.get_collection()) {
    StructureId our_sid = our_collection.structure_id;
    std::vector<ResidueIndex> filtered_indices;
    filtered_indices.reserve(our_collection.size());

    // for each residue in our structure
    for (ResidueIndex our_idx : our_collection.indices) {
      if (predicate(our_sid, our_idx)) { // apply th predicate to determine if this residue should be included
        filtered_indices.push_back(our_idx);
      }
    }

    // Add this structure's filtered residues to the result
    if (!filtered_indices.empty()) {
      result_collections.emplace_back(our_sid, std::move(filtered_indices));
    }
  }

  return ResidueSet(first.getMapping(), std::move(result_collections));
}
} // namespace detail


ResidueSet ResidueSet::operator|(const ResidueSet &other) const {
    detail::verify_mappings(mapping_, other.mapping_, "union");

    // create a lookup map: StructureId to indices
    std::unordered_map<StructureId, std::vector<ResidueIndex>> result_map;
    result_map.reserve(residue_collections_.size() + other.residue_collections_.size());

    // Seed with all collections from *this*
    for (auto const &collection : residue_collections_) {
        result_map.emplace(collection.structure_id, collection.indices);
    }

    // Merge in all collections from *other*
    for (auto const &other_collection : other.residue_collections_) {
        auto [it, inserted] = result_map.try_emplace(
            other_collection.structure_id,
            other_collection.indices  // only used if inserted == true
        );

        if (inserted) continue;

        auto &indices = it->second;
        std::vector<ResidueIndex> merged;
        merged.reserve(indices.size() + other_collection.indices.size());
        std::set_union(
            indices.begin(), indices.end(),
            other_collection.indices.begin(), other_collection.indices.end(),
            std::back_inserter(merged)
        );
        indices = std::move(merged);
    }

    std::vector<StructureResidues> result_collections;
    result_collections.reserve(result_map.size());
    for (auto &[sid, indices] : result_map) {
        result_collections.emplace_back(sid, std::move(indices));
    }

    return ResidueSet(mapping_, std::move(result_collections));
}

ResidueSet ResidueSet::operator&(const ResidueSet &other) const {
  detail::verify_mappings(mapping_, other.mapping_, "intersection");

  // Pre-compute all mapping matrices we'll need in one batch
  MappingMatrixBatchProvider provider(mapping_);
  auto matrices = provider.precompute_matrices(
      detail::extract_structure_ids(residue_collections_),
      detail::extract_structure_ids(other.residue_collections_)
  );

  return detail::rs_set_operation(*this, other, [this, &other, &matrices](StructureId our_sid, ResidueIndex our_idx) {
    // Check if *this* residue maps to any residue in the *other* set
    for (const auto &other_collection : other.residue_collections_) {
      StructureId other_sid = other_collection.structure_id;

      // same structure case
      if (detail::exists_in_same_structure(other_collection, our_sid, our_idx)) {
        return true;
      }

      // check mapping to *other* structures using pre-computed matrix
      auto matrix = matrices.get(our_sid, other_sid);
      if (detail::maps_to_any_in_collection(matrix, our_idx, other_collection)) {
        return true;
      }
    }
    return false;
  });
}

ResidueSet ResidueSet::operator-(const ResidueSet &other) const {
  detail::verify_mappings(mapping_, other.mapping_, "difference");

  // Pre-compute all mapping matrices we'll need in one batch
  MappingMatrixBatchProvider provider(mapping_);
  auto matrices = provider.precompute_matrices(
      detail::extract_structure_ids(residue_collections_),
      detail::extract_structure_ids(other.residue_collections_)
  );

  return detail::rs_set_operation(*this, other, [this, &other, &matrices](StructureId our_sid, ResidueIndex our_idx) {
    // Check if *this* residue maps to any residue in the *other* set
    for (const auto &other_collection : other.residue_collections_) {
      StructureId other_sid = other_collection.structure_id;

      if (detail::exists_in_same_structure(other_collection, our_sid, our_idx)) {
        return false;
      }

      auto matrix = matrices.get(our_sid, other_sid);
      if (detail::maps_to_any_in_collection(matrix, our_idx, other_collection)) {
        return false;
      }
    }
    return true;
  });
}

ResidueSet ResidueSet::intersection_without_batching(const ResidueSet &other) const {
  return detail::rs_set_operation(*this, other, [this, &other](StructureId our_sid, ResidueIndex our_idx) {
    // Check if *this* residue maps to any residue in the *other* set
    for (const auto &other_collection : other.residue_collections_) {

      if (detail::exists_in_same_structure(other_collection, our_sid, our_idx)) {
        return true;
      }

      if (detail::maps_to_any_in_collection(mapping_, our_sid, our_idx, other_collection)) {
        return true;
      }
    }
    return false;
  });
}

ResidueSet ResidueSet::difference_without_batching(const ResidueSet &other) const {
  return detail::rs_set_operation(*this, other, [this, &other](StructureId our_sid, ResidueIndex our_idx) {
    // Check if *this* residue maps to any residue in the *other* set
    for (const auto &other_collection : other.residue_collections_) {

      if (detail::exists_in_same_structure(other_collection, our_sid, our_idx)) {
        return false;
      }

      if (detail::maps_to_any_in_collection(mapping_, our_sid, our_idx, other_collection)) {
        return false;
      }
    }
    return true;
  });
}


namespace ResidueSetFactory {

ResidueSet empty(std::shared_ptr<BitSetMapping> mapping) { return ResidueSet(std::move(mapping), {}); }

ResidueSet from_structure(std::shared_ptr<BitSetMapping> mapping, StructureId structure_id) {

  // residue indices for *this* structure
  std::vector<ResidueIndex> indices = mapping->get_residue_indices(structure_id);

  // ResidueSet with these residues
  std::vector<ResidueSet::StructureResidues> collections;
  if (!indices.empty()) {
    collections.emplace_back(structure_id, std::move(indices));
  }

  return ResidueSet(std::move(mapping), std::move(collections));
}

ResidueSet from_structures(std::shared_ptr<BitSetMapping> mapping, const std::vector<StructureId> &structure_ids) {

  std::vector<ResidueSet::StructureResidues> collections;
  collections.reserve(structure_ids.size());

  // residue indices for each structure
  for (StructureId id : structure_ids) {
    std::vector<ResidueIndex> indices = mapping->get_residue_indices(id);
    if (!indices.empty()) {
      collections.emplace_back(id, std::move(indices));
    }
  }

  return ResidueSet(std::move(mapping), std::move(collections));
}

} // namespace ResidueSetFactory

} // namespace lahuta::mapping
