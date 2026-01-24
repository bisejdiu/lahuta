#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/GraphMol/Rings.h>

#include "aromatics.hpp"
#include "chemistry/common.hpp"
#include "logging/logging.hpp"
#include "planarity.hpp"
#include "residues/definitions.hpp"
#include "selections/mol_filters.hpp"

namespace lahuta {

// TODO: 1. since we use table-lookup to define (ideally most) aromatic residues, we should
//          also valide that they conform to expected aromatic ring properties
void add_rings_to_mol(const RDKit::RWMol &mol, const RDKit::VECT_INT_VECT &rings) {
  RDKit::VECT_INT_VECT bidx;
  RingUtils::convertToBonds(rings, bidx, mol);
  mol.getRingInfo()->addAllRings(rings, bidx);
}

AromaticRing get_molops_aromatic_rings(RDKit::RWMol &mol) {
  RDKit::VECT_INT_VECT rings, bonds;
  RDKit::MolOps::symmetrizeSSSR(mol, true);
  for (const std::vector<int> &ring : mol.getRingInfo()->atomRings()) {
    if (common::is_ring_aromatic(mol, ring)) {
      rings.push_back(ring);
    }

    // planarity check adapted from molstar: src/mol-model/structure/structure/unit/rings.ts (b67eda7)
    if (ring.size() < 5) continue;
    // avoid planarity check if ring contains any aromatic atoms
    if (common::has_any_aromatic_atom(mol, ring)) continue;
    if (!is_planar(mol, ring)) continue;
    rings.push_back(ring);
  }

  RingUtils::convertToBonds(rings, bonds, mol);
  return {rings, bonds};
}

std::vector<std::vector<int>>
map_rings(const std::vector<std::vector<int>> &aromatic_rings, const std::vector<int> &indices) {
  std::vector<std::vector<int>> mapped_rings;
  mapped_rings.reserve(aromatic_rings.size());

  for (const auto &ring : aromatic_rings) {
    std::vector<int> mapped_ring;
    mapped_ring.reserve(ring.size());
    for (int atom_idx : ring) {
      mapped_ring.push_back(indices[atom_idx]);
    }
    mapped_rings.push_back(std::move(mapped_ring));
  }

  return mapped_rings;
}

void apply_sssr_and_planarity_aromaticity(const RDKit::RWMol &mol, std::vector<int> &indices) {

  auto new_mol = filter_with_bonds(mol, indices);
  auto aromatic_rings = get_molops_aromatic_rings(new_mol);

  auto mapped_rings = map_rings(aromatic_rings.rings, indices);
  add_rings_to_mol(mol, mapped_rings);
}

void initialize_and_populate_ringinfo(const RDKit::RWMol &mol, const Residues &residues) {

  if (mol.getRingInfo()->isInitialized()) {
    mol.getRingInfo()->reset();
  }
  mol.getRingInfo()->initialize(RDKit::FIND_RING_TYPE_SYMM_SSSR);

  auto rings = find_and_process_aromatic_residues(mol, residues);
  add_rings_to_mol(mol, rings);

  // if (spdlog::should_log(spdlog::level::debug)) {
  if (true) { // FIX:
    auto unk_res = residues.filter(std::not_fn(definitions::is_predefined));

    std::unordered_map<std::string, int> residue_counts;
    for (const auto &res : unk_res) {
      residue_counts[res.name]++;
    }
    for (const auto &[name, count] : residue_counts) {
      Logger::get_logger()->debug("unk residue: {}, {}", name, count);
    }
  }

  auto unk_indices = residues.filter(std::not_fn(definitions::is_predefined)).get_atom_ids();
  if (!unk_indices.empty()) {
    std::sort(unk_indices.begin(), unk_indices.end());
    apply_sssr_and_planarity_aromaticity(mol, unk_indices);
  }
}

} // namespace lahuta
