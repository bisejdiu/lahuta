#include "contacts/contacts.hpp"
#include "lahuta.hpp"

namespace lahuta {

std::vector<const RDKit::Atom *> get_atom_types(const Luni *luni, AtomType type) {
  const auto &atom_types = luni->get_atom_types();
  const auto &mol = luni->get_molecule();
  std::vector<const RDKit::Atom *> filtered_atoms;
  for (const auto *atom : mol.atoms()) {
    /*if (AtomTypeFlags::has(atom_types[atom->getIdx()], type)) {*/
    if (AtomTypeFlags::has_any(atom_types[atom->getIdx()], type)) {
      filtered_atoms.push_back(atom);
    }
  }
  return filtered_atoms;
}

// FIX: this could be a Luna method
const AtomDataVec
get_atom_data(const Luni *luni, AtomType type, FeatureTypeCheckFunc check_func) {
  AtomDataVec atom_data_vec;
  const auto &atom_types = luni->get_atom_types();
  const auto &mol = luni->get_molecule();
  const auto &conf = mol.getConformer();

  for (const auto *atom : mol.atoms()) {
    /*if (AtomTypeFlags::has_any(atom_types[atom->getIdx()], type)) {*/
    if (check_func(atom_types[atom->getIdx()], type)) {
      /*atom_data_vec.data.emplace_back(*/
      /*    atom_types[atom->getIdx()], atom, &conf.getAtomPos(atom->getIdx()), atom->getIdx());*/
      atom_data_vec.add_data(mol, atom, atom_types[atom->getIdx()]);
    }
  }
  return atom_data_vec;
}

} // namespace lahuta
