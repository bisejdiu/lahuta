#include "contacts/aromaticity.hpp"
#include "typing/types.hpp"

namespace lahuta {

// FIX: this may not be needed any longer!
std::vector<const RDKit::Atom *> get_atoms(const RDKit::RWMol &mol, const std::vector<int> &indices) {
  std::vector<const RDKit::Atom *> atoms;
  atoms.reserve(indices.size());
  for (int idx : indices) {
    atoms.push_back(mol.getAtomWithIdx(idx));
  }
  return atoms;
}

std::vector<GroupRec> add_aromatic_rings(const RDKit::RWMol &mol, const Residues &) {
  std::vector<GroupRec> groups;

  for (const auto &ring : mol.getRingInfo()->atomRings()) {
    auto is_aromatic = [&mol](int idx) { return mol.getAtomWithIdx(idx)->getIsAromatic(); };
    if (std::all_of(ring.begin(), ring.end(), is_aromatic)) {

      std::vector<std::reference_wrapper<const RDKit::Atom>> atoms;
      atoms.reserve(ring.size());

      for (int idx : ring) {
        atoms.push_back(std::cref(*mol.getAtomWithIdx(idx)));
      }

      groups.push_back(GroupRec{
        /*.a_type =*/ AtomType::Aromatic,
        /*.type   =*/ FeatureGroup::None,
        /*.atoms  =*/ std::move(atoms),
        /*.center =*/ RDGeom::Point3D(0,0,0)
      });
    }
  }

  return groups;
}

} // namespace lahuta
