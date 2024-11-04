#include "contacts/hydrophobic.hpp"
#include "contacts/search.hpp"
#include "entities.hpp"
#include "lahuta.hpp"

namespace lahuta {

AtomType add_hydrophobic_atom(const RDKit::RWMol &mol, const RDKit::Atom &atom) {
  const int atomic_number = atom.getAtomicNum();

  if (atomic_number == 9) return AtomType::HYDROPHOBIC;
  if (atomic_number != 6) return AtomType::NONE;

  for (const auto &bond : mol.atomBonds(&atom)) {
    const RDKit::Atom *neighbor = bond->getOtherAtom(&atom);
    int neighbor_atomic_num = neighbor->getAtomicNum();

    if (neighbor_atomic_num != 6 && neighbor_atomic_num != 1) {
      return AtomType::NONE;
    }
  }

  return AtomType::HYDROPHOBIC;
}

Contacts find_hydrophobic_bonds(const Luni &luni, HydrophobicParams opts) {

  Contacts contacts(&luni);
  const auto hydrophobic_atoms = AtomEntityCollection::filter(&luni, AtomType::HYDROPHOBIC);

  EntityNeighborSearch ens(luni.get_conformer());
  NSResults results = ens.search(hydrophobic_atoms, opts.distance_max);

  for (const auto &[pair, dist] : results) {
    auto [atom1_index, atom2_index] = pair;
    const auto &atom1_data = hydrophobic_atoms.get_data()[atom1_index];
    const auto &atom2_data = hydrophobic_atoms.get_data()[atom2_index];

    if (are_residueids_close(luni.get_molecule(), *atom1_data.atom, *atom2_data.atom, 0)) continue;
    if (atom1_data.atom->getAtomicNum() == 9 && atom2_data.atom->getAtomicNum() == 9) continue;

    contacts.add(Contact(
        static_cast<EntityID>(atom1_data.atom->getIdx()),
        static_cast<EntityID>(atom2_data.atom->getIdx()),
        dist,
        InteractionType::Hydrophobic));
  }

  return contacts;
}

} // namespace lahuta
