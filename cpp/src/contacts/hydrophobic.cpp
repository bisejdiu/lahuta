#include "contacts/hydrophobic.hpp"
#include "contacts/search.hpp"
#include "contacts/utils.hpp"
#include "entities.hpp"
#include "lahuta.hpp"
#include "elements.hpp"

namespace lahuta {

AtomType add_hydrophobic_atom(const RDKit::RWMol &mol, const RDKit::Atom &atom) {
  const auto at_n = atom.getAtomicNum();

  if (at_n == Element::F) return AtomType::HYDROPHOBIC;
  if (at_n != Element::C) return AtomType::NONE;

  for (const auto &bond : mol.atomBonds(&atom)) {
    const RDKit::Atom *nbr = bond->getOtherAtom(&atom);
    const auto nbr_at_n = nbr->getAtomicNum();
    if (nbr_at_n != Element::C && nbr_at_n != Element::H) return AtomType::NONE;
  }

  return AtomType::HYDROPHOBIC;
}

Contacts find_hydrophobic_bonds(const Luni &luni, std::optional<HydrophobicParams> params) {

  Contacts contacts(&luni);
  HydrophobicParams opts = params.value_or(HydrophobicParams{});

  const auto hydrophobic_atoms = AtomEntityCollection::filter(&luni, AtomType::HYDROPHOBIC);
  NSResults results = EntityNeighborSearch::search(hydrophobic_atoms, opts.distance_max);

  for (const auto &[pair, dist] : results) {
    auto [atom1_index, atom2_index] = pair;
    const auto &atom1_data = hydrophobic_atoms.get_data()[atom1_index];
    const auto &atom2_data = hydrophobic_atoms.get_data()[atom2_index];

    if (are_residueids_close(luni.get_molecule(), *atom1_data.atom, *atom2_data.atom, 0)) continue;
    if (atom1_data.atom->getAtomicNum() == Element::F && atom2_data.atom->getAtomicNum() == Element::F) continue;

    contacts.add(Contact(
        static_cast<EntityID>(atom1_data.atom->getIdx()),
        static_cast<EntityID>(atom2_data.atom->getIdx()),
        dist,
        InteractionType::Hydrophobic));
  }

  return contacts;
}

} // namespace lahuta
