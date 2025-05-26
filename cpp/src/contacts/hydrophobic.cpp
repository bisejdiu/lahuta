#include "contacts/hydrophobic.hpp"
#include "elements.hpp"
#include "typing/types.hpp"

namespace lahuta {

AtomType add_hydrophobic_atom(const RDKit::RWMol &mol, const RDKit::Atom &atom) {
  const auto at_n = atom.getAtomicNum();

  if (at_n == Element::F) return AtomType::Hydrophobic;
  if (at_n != Element::C) return AtomType::None;

  for (const auto &bond : mol.atomBonds(&atom)) {
    const RDKit::Atom *nbr = bond->getOtherAtom(&atom);
    const auto nbr_at_n = nbr->getAtomicNum();
    if (nbr_at_n != Element::C && nbr_at_n != Element::H) return AtomType::None;
  }

  return AtomType::Hydrophobic;
}

// ContactSet find_hydrophobic_bonds(const Topology& topology, const HydrophobicParams& opts) {
//   return find_contacts(
//     {topology, opts},
//     [](const AtomRec &rec) { return (rec.type & AtomType::Hydrophobic) == AtomType::Hydrophobic; },
//     {opts.distance_max, 0.4},
//     [](std::uint32_t rec_idx_a, std::uint32_t rec_idx_b, float dist, const ContactContext& ctx) -> InteractionType {
//       const auto& opts = ctx.get_params<HydrophobicParams>();
//       const auto &mol = ctx.topology.molecule();
//
//       const auto &atom_a = ctx.topology.atom(rec_idx_a).atom;
//       const auto &atom_b = ctx.topology.atom(rec_idx_b).atom;
//
//       if (are_residueids_close(mol, atom_a, atom_b, 0)) return InteractionType::None;
//       if (atom_a.getAtomicNum() == Element::F && atom_b.getAtomicNum() == Element::F) return InteractionType::None;
//
//       return InteractionType::Hydrophobic;
//     }
//   );
// }

} // namespace lahuta
