#ifndef LAHUTA_MAPPING_RESIDUE_SET_HPP
#define LAHUTA_MAPPING_RESIDUE_SET_HPP

#include "_defs.hpp"
#include <algorithm>
#include <memory>
#include <type_traits>
#include <vector>

// clang-format off
namespace lahuta::mapping {
class BitSetMapping;

// Provides set-like ops for collections of residues across multiple structures.
// Set operations work across structures by default.
class ResidueSet {
public:
  /// Stores a collection of residues from a single structure
  struct StructureResidues {
    StructureId structure_id;
    std::vector<ResidueIndex> indices; // sorted

    explicit StructureResidues(StructureId structure_id_, std::vector<ResidueIndex> indices_)
        : structure_id(structure_id_), indices(std::move(indices_)) {}

    /// number of residues in this structure
    size_t size() const { return indices.size(); }

    /// Check if this structure contains a specific residue index
    bool contains(ResidueIndex index) const { return std::binary_search(indices.begin(), indices.end(), index); }
  };

private:
  std::shared_ptr<BitSetMapping> mapping_;
  std::vector<StructureResidues> residue_collections_;
  size_t total_count_{0}; // Track the total number of residues

public:
  /// Create a ResidueSet with the specified mapping and residue collections
  explicit ResidueSet(std::shared_ptr<BitSetMapping> mapping, std::vector<StructureResidues> residue_collections)
      : mapping_(std::move(mapping)), residue_collections_(std::move(residue_collections)) {
    for (const auto &collection : residue_collections_) {
      total_count_ += collection.size();
    }
  }

  /// Get all structure-residue collections
  const std::vector<StructureResidues> &get_collection() const { return residue_collections_; }

  /// Get the total number of residues across all structures
  size_t size() const { return total_count_; }

  /// Check if the set is empty
  bool empty() const { return residue_collections_.empty() || total_count_ == 0; }

  /// Union:        A | B: residues are pooled from all structures, regardless of mapping
  /// Intersection: A & B: residues A that have equivalent residues in B
  /// Difference:   A - B: residues A that don't have equivalent residues in B
  /// Sym. diff.:   A ^ B: residues A or B that don't have equivalent residues in the other
  ResidueSet operator|(const ResidueSet &other) const;
  ResidueSet operator&(const ResidueSet &other) const;
  ResidueSet operator-(const ResidueSet &other) const;
  ResidueSet operator^(const ResidueSet &other) const {
    return (*this - other) | (other - *this);
  }

  ResidueSet &operator|=(const ResidueSet &other) { *this = *this | other; return *this; }
  ResidueSet &operator&=(const ResidueSet &other) { *this = *this & other; return *this; }
  ResidueSet &operator-=(const ResidueSet &other) { *this = *this - other; return *this; }
  ResidueSet &operator^=(const ResidueSet &other) { *this = *this ^ other; return *this; }

  /// Intersection without matrix batching (slower performance)
  ResidueSet intersection_without_batching(const ResidueSet &other) const;

  /// Difference without matrix batching (slower performance)
  ResidueSet difference_without_batching(const ResidueSet &other) const;

  /// Apply a function to each residue in the set
  template <typename Func>
  void for_each(Func &&func) const {
    static_assert(std::is_invocable_r_v<void, Func, StructureId, ResidueIndex>,
                  "Func must be callable as void(StructureId,ResidueIndex)");

    for (const auto &collection : residue_collections_) {
      StructureId sid = collection.structure_id;
      for (ResidueIndex index : collection.indices) {
        func(sid, index);
      }
    }
  }

  // Map each residue to a value using a mapping function and return a vector of results
  template <typename T, typename MapFunc>
  std::vector<T> map(MapFunc &&map_func) const {
    static_assert(std::is_invocable_r_v<T, MapFunc, StructureId, ResidueIndex>,
                  "MapFunc must be callable as T(StructureId,ResidueIndex)");

    std::vector<T> result; result.reserve(total_count_);
    for_each([&](StructureId sid, ResidueIndex index) { result.push_back(map_func(sid, index)); });

    return result;
  }

  // Filter residues based on a predicate function
  template <typename Predicate>
  ResidueSet filter(Predicate &&predicate) const {
    static_assert(std::is_invocable_r_v<bool, Predicate, StructureId, ResidueIndex>,
                  "Predicate must be callable as bool(StructureId,ResidueIndex)");

    std::vector<StructureResidues> filtered_collections;

    for (const auto &collection : residue_collections_) {
      StructureId sid = collection.structure_id;
      std::vector<ResidueIndex> filtered_indices;

      for (ResidueIndex index : collection.indices) {
        if (predicate(sid, index)) {
          filtered_indices.push_back(index);
        }
      }

      if (!filtered_indices.empty()) {
        filtered_collections.emplace_back(sid, std::move(filtered_indices));
      }
    }

    return ResidueSet(mapping_, std::move(filtered_collections));
  }

  /// Get the underlying mapping
  std::shared_ptr<BitSetMapping> getMapping() const { return mapping_; }
};

namespace ResidueSetFactory {

  /// Empty ResidueSet with just a mapping
  ResidueSet empty(std::shared_ptr<BitSetMapping> mapping);

  /// With all residues from a single structure
  ResidueSet from_structure(std::shared_ptr<BitSetMapping> mapping, StructureId structure_id);

  /// With all residues from a list of structures
  ResidueSet from_structures(std::shared_ptr<BitSetMapping> mapping, const std::vector<StructureId> &structure_ids);
}

} // namespace lahuta::mapping

#endif // LAHUTA_MAPPING_RESIDUE_SET_HPP
