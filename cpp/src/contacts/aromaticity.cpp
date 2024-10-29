#include "contacts/aromaticity.hpp"
#include "atom_types.hpp"

namespace lahuta {

std::vector<const RDKit::Atom *> get_atoms(const RDKit::RWMol &mol, const std::vector<int> &indices) {
  std::vector<const RDKit::Atom *> atoms;
  for (int idx : indices) {
    atoms.push_back(mol.getAtomWithIdx(idx));
  }
  return atoms;
}

FeatureVec add_aromatic_rings(const RDKit::RWMol &mol, const Residues &residues) {
  FeatureVec features;

  // FIX: here we check if RingInfo is initialized ourselves, and we should decide to either throw or run whatever ring computation logic we end up implementing

  for (const auto &ring : mol.getRingInfo()->atomRings()) {
    if (std::all_of(ring.begin(), ring.end(), [&mol](int idx) {
          return mol.getAtomWithIdx(idx)->getIsAromatic();
        })) {
      auto atoms = get_atoms(mol, ring);
      features.add_feature(create_feature(AtomType::AROMATIC, FeatureGroup::None, atoms));
    }
  }

  return features;
}

} // namespace lahuta
