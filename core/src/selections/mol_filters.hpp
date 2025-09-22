#pragma once

#include "logging.hpp"
#include <gemmi/mmread_gz.hpp>
#include <rdkit/Geometry/point.h>
#include <rdkit/GraphMol/MonomerInfo.h>
#include <rdkit/GraphMol/RWMol.h>
#include <vector>

namespace lahuta {

inline RDKit::RWMol filter_with_conf(RDKit::RWMol &mol, std::vector<int> &indices) {
  RDKit::Conformer conf = mol.getConformer();
  RDKit::RWMol new_mol;
  RDKit::Conformer *new_conf = new RDKit::Conformer();

  for (auto atomIdx : indices) {
    auto atom = mol.getAtomWithIdx(atomIdx);
    RDKit::Atom *newAtom = new RDKit::Atom(*atom);
    new_mol.addAtom(newAtom, true, true);
    auto pos = conf.getAtomPos(atom->getIdx());
    new_conf->setAtomPos(newAtom->getIdx(), pos);
  }
  new_mol.addConformer(new_conf, true);
  // new_mol.updatePropertyCache(false);

  return new_mol;
}

// indices will be sorted in ascending order
inline RDKit::RWMol filter_with_bonds(const RDKit::RWMol &mol, std::vector<int> &indices) {

  RDKit::RWMol filtered_mol;
  const RDKit::Conformer &conf = mol.getConformer();

  auto *filtered_conf = new RDKit::Conformer();

  std::unordered_map<int, int> atom_index_map;

  if (indices.empty()) {
    filtered_mol.addConformer(filtered_conf, true);
    filtered_mol.updatePropertyCache(false);
    return filtered_mol;
  }

  std::sort(indices.begin(), indices.end());

  if (indices.back() >= mol.getNumAtoms()) {
    throw std::out_of_range("Index out of range");
  }

  // copy atoms and positions.
  for (size_t i = 0; i < indices.size(); ++i) {
    int atom_idx = indices[i];
    const RDKit::Atom *atom = mol.getAtomWithIdx(atom_idx);
    int filtered_idx = filtered_mol.addAtom(new RDKit::Atom(*atom), false, true);
    atom_index_map[atom_idx] = filtered_idx;
    filtered_conf->setAtomPos(filtered_idx, conf.getAtomPos(atom_idx));
  }

  // add bonds only if both atoms are in the selected indices.
  for (const auto &bond : mol.bonds()) {
    int start_idx = bond->getBeginAtomIdx();
    int end_idx = bond->getEndAtomIdx();

    if (atom_index_map.find(start_idx) != atom_index_map.end() && atom_index_map.find(end_idx) != atom_index_map.end()) {
      int new_start_idx = atom_index_map[start_idx];
      int new_end_idx = atom_index_map[end_idx];
      filtered_mol.addBond(new_start_idx, new_end_idx, bond->getBondType());
    }
  }

  filtered_mol.addConformer(filtered_conf, true);
  filtered_mol.updatePropertyCache(false);

  Logger::get_logger()->debug(
      "Filtered molecule has {} atoms and {} bonds",
      filtered_mol.getNumAtoms(),
      filtered_mol.getNumBonds());

  return filtered_mol;
}

} // namespace lahuta
