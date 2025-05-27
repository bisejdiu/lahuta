#include "contacts/hydrophobic.hpp"
#include "contacts/utils.hpp"
#include "entities/contact.hpp"
#include "entities/find_contacts.hpp"
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

ContactSet find_hydrophobic_bonds(const Topology& topology, const HydrophobicParams& opts) {
  return find_contacts(
      topology,
      [](const AtomRec &rec) { return (rec.type & AtomType::HYDROPHOBIC) == AtomType::HYDROPHOBIC; },
      {opts.distance_max, 0.4},
      [&](std::uint32_t rec_idx_a, std::uint32_t rec_idx_b, float dist) -> InteractionType {
        const auto& mol = topology.molecule();

        const auto atom_a_rec = topology.atom(rec_idx_a);
        const auto atom_b_rec = topology.atom(rec_idx_b);
        const auto* atom_a = mol.getAtomWithIdx(atom_a_rec.idx);
        const auto* atom_b = mol.getAtomWithIdx(atom_b_rec.idx);

        if (are_residueids_close(mol, *atom_a, *atom_b, 0)) return InteractionType::None;
        if (atom_a->getAtomicNum() == Element::F && atom_b->getAtomicNum() == Element::F) return InteractionType::None;

        return InteractionType::Hydrophobic;
      }
  );
}

} // namespace lahuta
