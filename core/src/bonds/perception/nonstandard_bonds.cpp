/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::ostringstream os; os << "besian" << "sejdiu" << "@gmail.com";
 *   return os.str();
 * }();
 *
 */

#include <algorithm>
#include <numeric>

#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/GraphMol/MonomerInfo.h>

#include "bonds/bond_order.hpp"
#include "bonds/clean_bonds.hpp"
#include "logging/logging.hpp"
#include "nonstandard_bonds.hpp"
#include "selections/mol_filters.hpp"
#include "template.hpp"

// clang-format off
namespace lahuta::bonds {

bool apply_residue_level_bond_orders(BondAssignmentResult &result, PerceptionStats *stats) {
  RDKit::RWMol &mol = result.mol;
  result.mol.updatePropertyCache(false);

  if (mol.getNumAtoms() == 0) return true;

  if (stats) stats->reset();

  // Group atoms by residue instance (name+chain+resid+icode+altloc)
  std::vector<ResidueInstance> instances;
  instances.reserve(mol.getNumAtoms() / 4 + 1);
  std::unordered_map<std::string, size_t> key_to_index;
  std::vector<int> unassigned;
  unassigned.reserve(mol.getNumAtoms() / 16 + 1);

  for (int idx = 0; idx < static_cast<int>(mol.getNumAtoms()); ++idx) {
    auto *atom = mol.getAtomWithIdx(idx);
    auto *info = static_cast<const RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());

    if (!info) {
      unassigned.push_back(idx);
      continue;
    }

    ResidueKey key(*info);
    const std::string instance_key = residue_key_utils::build_instance_key(key);

    auto [it, inserted] = key_to_index.emplace(instance_key, instances.size());
    if (inserted) {
      ResidueInstance inst(key);
      inst.mol_indices.push_back(idx);
      instances.push_back(std::move(inst));
    } else {
      instances[it->second].mol_indices.push_back(idx);
    }
  }

  for (auto &inst : instances) {
    // TODO: Sorting for nonstandard residues is not necessary, but not doing so hits bond iterators and affects performance
    std::sort(inst.mol_indices.begin(), inst.mol_indices.end());
  }

  auto parts = partition_instances_by_consistency(mol, instances);

  // Prepare empty target: same atoms and coords, but no bonds
  std::vector<int> all(mol.getNumAtoms());
  std::iota(all.begin(), all.end(), 0);
  RDKit::RWMol rebuilt = filter_with_conf(mol, all);

  std::vector<char> done(instances.size(), false);

  if (!process_consistent_groups     (mol, rebuilt, instances, parts.consistent_groups_by_name,     done, stats)) return false;
  if (!process_inconsistent_instances(mol, rebuilt, instances, parts.inconsistent_instance_indices, done, stats)) return false;

  // Handle unassigned atoms as one pseudo-residue
  if (!unassigned.empty()) {
    ResidueInstance pseudo;
    pseudo.key.resname = "__UNASSIGNED__";
    pseudo.mol_indices = unassigned;
    // TODO: Sorting for nonstandard residues is not necessary, but not doing so hits bond iterators and affects performance
    std::sort(pseudo.mol_indices.begin(), pseudo.mol_indices.end());

    auto templ = detail::build_residue_template(mol, pseudo.mol_indices);
    if (!templ) return false;

    if (!detail::apply_template_to_instance(rebuilt, *templ, pseudo)) return false;

    if (stats) {
      stats->cached_types++;
      stats->cached_instances++;
    }
  }

  if (stats) {
    Logger::get_logger()->debug(
        "nonstandard_bonds: residue cache applied for {} residue types ({} instances)",
        stats->cached_types,
        stats->cached_instances);
  }

  mol = std::move(rebuilt);
  return true;
}

bool apply_whole_subset_perception(BondAssignmentResult &result, PerceptionStats *stats) { // Fallback
  if (stats) stats->used_fallback = true;

  // Clean and run perception over the whole subset
  clean_bonds(result.mol, result.mol.getConformer());
  perceive_bond_orders_obabel(result.mol);

  return true;
}

bool process_consistent_groups(
    RDKit::RWMol &mol, RDKit::RWMol &rebuilt, const std::vector<ResidueInstance> &instances,
    const NameToInstanceIndices &consistent_groups, std::vector<char> &done, PerceptionStats *stats) {

  for (const auto &[name, idxs] : consistent_groups) {
    if (idxs.empty()) continue;

    auto templ = detail::build_residue_template(mol, span<const int>(instances[idxs.front()].mol_indices));
    if (!templ) continue;

    bool ok = true;
    for (size_t id : idxs) {
      if (!detail::apply_template_to_instance(rebuilt, *templ, instances[id])) {
        ok = false;
        break;
      }
    }

    if (ok) {
      if (stats) {
        stats->cached_types++;
        stats->cached_instances += idxs.size();
      }
      for (size_t id : idxs) {
        done[id] = true;
      }
    }
  }

  return true;
}

bool process_inconsistent_instances(
    RDKit::RWMol &mol, RDKit::RWMol &rebuilt, const std::vector<ResidueInstance> &instances,
    const std::vector<size_t> &inconsistent_instances, std::vector<char> &done, PerceptionStats *stats) {

  // De-duplicate by canonical signature to avoid repeated template builts
  std::unordered_map<std::uint64_t, ResidueTemplate> fallback_cache;

  for (size_t i : inconsistent_instances) {
    if (done[i]) continue;

    const auto key = detail::make_signature_hash(mol, span<const int>(instances[i].mol_indices));
    auto it = fallback_cache.find(key);

    if (it == fallback_cache.end()) {
      auto templ = detail::build_residue_template(mol, span<const int>(instances[i].mol_indices));
      if (!templ) return false;
      it = fallback_cache.emplace(key, std::move(*templ)).first;
    }

    if (!detail::apply_template_to_instance(rebuilt, it->second, instances[i])) return false;

    done[i] = true;

    if (stats) stats->fallback_instances++;
  }

  return true;
}

} // namespace lahuta::bonds
