#include "contacts/aromaticity.hpp"
#include "GraphMol/Rings.h"
#include "RDGeneral/types.h"
#include "atom_types.hpp"
#include "contacts/charges.hpp"
#include "convert.hpp"

namespace lahuta {

void add_rings_to_mol(const RDKit::RWMol &mol, const RDKit::VECT_INT_VECT &rings) {
  RDKit::VECT_INT_VECT bidx;
  RingUtils::convertToBonds(rings, bidx, mol);
  for (size_t i = 0; i < rings.size(); ++i) {
    mol.getRingInfo()->addRing(rings[i], bidx[i]);
  }
}

std::vector<const RDKit::Atom *> get_atoms(const RDKit::RWMol &mol, const std::vector<int> &indices) {
  std::vector<const RDKit::Atom *> atoms;
  for (int idx : indices) {
    atoms.push_back(mol.getAtomWithIdx(idx));
  }
  return atoms;
}

AromaticRing get_molops_aromatic_rings(RDKit::RWMol &mol) {
  RDKit::VECT_INT_VECT rings, bonds;
  RDKit::MolOps::symmetrizeSSSR(mol, true);
  for (const auto &ring : mol.getRingInfo()->atomRings()) {
    if (std::all_of(ring.begin(), ring.end(), [&mol](int idx) {
          return mol.getAtomWithIdx(idx)->getIsAromatic();
        })) {
      rings.push_back(ring);
    }
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

FeatureVec add_aromatic_rings(const RDKit::RWMol &mol, const std::vector<int> &indices) {
  FeatureVec features;

  // create new molecule with only the selected atoms
  auto new_mol = filter_with_bonds(mol, indices);
  auto aromatic_rings = get_molops_aromatic_rings(new_mol);

  auto mapped_rings = map_rings(aromatic_rings.rings, indices);
  auto mapped_bonds = map_rings(aromatic_rings.bonds, indices);

  for (size_t i = 0; i < mapped_rings.size(); ++i) {
    // populate ring info
    mol.getRingInfo()->addRing(mapped_rings[i], mapped_bonds[i]);

    // populate features
    auto atoms = get_atoms(mol, mapped_rings[i]);
    features.add_feature(create_feature(AtomType::AROMATIC, FeatureGroup::None, atoms));
  }

  return features;
}

FeatureVec add_aromatic_rings(const RDKit::RWMol &mol, Residues &residues) {
  FeatureVec features;

  auto rings = residue_props::find_aromatic_rings(mol, residues);
  add_rings_to_mol(mol, rings);

  for (const auto &ring : rings) {
    auto atoms = get_atoms(mol, ring);
    features.add_feature(create_feature(AtomType::AROMATIC, FeatureGroup::None, atoms));
  }

  // compute rings for unknown residues
  auto unk_indices = residue_props::get_unknown_residues<std::vector<int>>(residues, AminoAcidNames);
  if (!unk_indices.empty()) {
    std::sort(unk_indices.begin(), unk_indices.end());
    auto additional_features = add_aromatic_rings(mol, unk_indices);
    features.merge(additional_features);
  }

  return features;
}

} // namespace lahuta
