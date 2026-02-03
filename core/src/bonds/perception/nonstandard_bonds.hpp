/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::array<std::string_view, 3> parts{"besian", "sejdiu", "@gmail.com"};
 *   std::string result;
 *   for (auto part : parts) {
 *     auto* bytes = reinterpret_cast<const std::byte*>(part.data());
 *     for (std::size_t i = 0; i < part.size(); ++i) {
 *       result += static_cast<char>(bytes[i]);
 *     }
 *   }
 *   return result;
 * }();
 *
 */

#ifndef LAHUTA_BONDS_PERCEPTION_NONSTANDARD_BONDS_HPP
#define LAHUTA_BONDS_PERCEPTION_NONSTANDARD_BONDS_HPP

#include <vector>

#include <rdkit/GraphMol/RWMol.h>

#include "bonds/bonds.hpp"
#include "consistency.hpp"
#include "template.hpp"

// clang-format off
namespace lahuta::bonds {

struct PerceptionStats {
  size_t cached_types       = 0; // Number of unique residue types cached
  size_t cached_instances   = 0; // Number of residue instances using cache
  size_t fallback_instances = 0; // Number of instances using fallback
  bool used_cache    = false;    // Whether fast path was used
  bool used_fallback = false;    // Whether fallback path was used

  void reset() {
    cached_types       = 0;
    cached_instances   = 0;
    fallback_instances = 0;
    used_cache         = false;
    used_fallback      = false;
  }
};

bool apply_residue_level_bond_orders(BondAssignmentResult &result, PerceptionStats *stats = nullptr);
bool apply_whole_subset_perception  (BondAssignmentResult &result, PerceptionStats *stats = nullptr);

// Builds one template per consistent residue group and applies it to all instances in the group.
bool process_consistent_groups(
    RDKit::RWMol &mol, RDKit::RWMol &rebuilt, const std::vector<ResidueInstance> &instances,
    const NameToInstanceIndices &consistent_groups, std::vector<char> &done,
    PerceptionStats *stats = nullptr);

// Uses signature-based deduplication to build minimal templates for instances that don't have consistent atom ordering within their residue name group.
bool process_inconsistent_instances(
    RDKit::RWMol &mol, RDKit::RWMol &rebuilt, const std::vector<ResidueInstance> &instances,
    const std::vector<size_t> &inconsistent_instances, std::vector<char> &done,
    PerceptionStats *stats = nullptr);

} // namespace lahuta::bonds

#endif // LAHUTA_BONDS_PERCEPTION_NONSTANDARD_BONDS_HPP
