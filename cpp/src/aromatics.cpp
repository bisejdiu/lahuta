#include "aromatics.hpp"
#include "GraphMol/MolOps.h"
#include "GraphMol/Rings.h"
#include "common.hpp"
#include "convert.hpp"
#include "definitions.hpp"
#include "planarity.hpp"

namespace lahuta {

void add_rings_to_mol(const RDKit::RWMol &mol, const RDKit::VECT_INT_VECT &rings) {
  RDKit::VECT_INT_VECT bidx;
  RingUtils::convertToBonds(rings, bidx, mol);
  for (size_t i = 0; i < rings.size(); ++i) {
    mol.getRingInfo()->addRing(rings[i], bidx[i]);
  }
}

AromaticRing get_molops_aromatic_rings(RDKit::RWMol &mol) {
  RDKit::VECT_INT_VECT rings, bonds;
  RDKit::MolOps::symmetrizeSSSR(mol, true);
  for (const std::vector<int> &ring : mol.getRingInfo()->atomRings()) {
    if (common::is_ring_aromatic(mol, ring)) {
      rings.push_back(ring);
    }

    // planarity check adapted from
    // https://github.com/molstar/molstar/blob/master/src/mol-model/structure/structure/unit/rings.ts
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

void apply_sssr_and_planarity_aromaticity(const RDKit::RWMol &mol, const std::vector<int> &indices) {

  // create new molecule with only the selected atoms
  auto new_mol = filter_with_bonds(mol, indices);
  auto aromatic_rings = get_molops_aromatic_rings(new_mol);

  auto mapped_rings = map_rings(aromatic_rings.rings, indices);
  auto mapped_bonds = map_rings(aromatic_rings.bonds, indices);

  add_rings_to_mol(mol, mapped_rings);
}

void initialize_and_populate_ringinfo(const RDKit::RWMol &mol, const Residues &residues) {
  using namespace residue_props;

  if (mol.getRingInfo()->isInitialized()) {
    mol.getRingInfo()->reset();
  }
  mol.getRingInfo()->initialize(RDKit::FIND_RING_TYPE_SYMM_SSSR);

  auto rings = residue_props::tbl_find_aromatic_rings(mol, residues);
  add_rings_to_mol(mol, rings);

  // compute rings for unknown residues
  auto unk_indices = get_unknown_residues<std::vector<int>>(residues, definitions::is_protein_extended);
  if (!unk_indices.empty()) {
    std::sort(unk_indices.begin(), unk_indices.end());
    apply_sssr_and_planarity_aromaticity(mol, unk_indices);
  }
}

} // namespace lahuta
