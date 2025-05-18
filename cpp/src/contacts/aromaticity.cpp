#include "contacts/aromaticity.hpp"
#include "atom_types.hpp"

namespace lahuta {

std::vector<const RDKit::Atom *> get_atoms(const RDKit::RWMol &mol, const std::vector<int> &indices) {
  std::vector<const RDKit::Atom *> atoms;
  atoms.reserve(indices.size());
  for (int idx : indices) {
    atoms.push_back(mol.getAtomWithIdx(idx));
  }
  return atoms;
}

GroupEntityCollection add_aromatic_rings(const RDKit::RWMol &mol, const Residues &residues) {
  GroupEntityCollection features;

  for (const auto &ring : mol.getRingInfo()->atomRings()) {
    auto is_aromatic = [&mol](int idx) { return mol.getAtomWithIdx(idx)->getIsAromatic(); };
    if (std::all_of(ring.begin(), ring.end(), is_aromatic)) {
      auto atoms = get_atoms(mol, ring);
      features.add_data(AtomType::AROMATIC, FeatureGroup::None, atoms);
    }
  }

  return features;
}

} // namespace lahuta
