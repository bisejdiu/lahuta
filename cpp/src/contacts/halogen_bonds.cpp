#include "contacts/halogen_bonds.hpp"
/*#include "contacts/halo_geo_validity.hpp"*/
#include "entities/find_contacts.hpp"
#include "entities/context.hpp"
#include <elements.hpp>

namespace lahuta {

std::unordered_set<Element> HalogenDonors    = {Element::Cl, Element::Br, Element::I};
std::unordered_set<Element> HalogenAcceptors = {Element::N,  Element::O,  Element::S};
std::unordered_set<Element> HalogenBinders   = {Element::C,  Element::N,  Element::P, Element::S};

AtomType add_halogen_donor(const RDKit::RWMol &mol, const RDKit::Atom &atom) {
  const auto at_n = static_cast<Element>(atom.getAtomicNum());
  if (HalogenDonors.count(at_n)) return AtomType::XbondDonor;
  return AtomType::None;
}

AtomType add_halogen_acceptor(const RDKit::RWMol &mol, const RDKit::Atom &atom) {
  const auto at_n = static_cast<Element>(atom.getAtomicNum());
  if (!HalogenAcceptors.count(at_n)) return AtomType::None;

  for (const auto bond : mol.atomBonds(&atom)) {
    const RDKit::Atom *nbr = bond->getOtherAtom(&atom);
    const auto nbr_at_n = static_cast<Element>(nbr->getAtomicNum());

    if (HalogenBinders.count(nbr_at_n)) return AtomType::XBondAcceptor;
  }

  return AtomType::None;
}

// ContactSet find_halogen_bonds(const Topology &topology, const HalogenParams &params) {
//   return find_contacts(
//     {topology, params},
//     [](const AtomRec &rec) { return (rec.type & AtomType::XbondDonor)    == AtomType::XbondDonor; },
//     [](const AtomRec &rec) { return (rec.type & AtomType::XBondAcceptor) == AtomType::XBondAcceptor; },
//     {params.distance_max, 0, 0, 0.7},
//     [](std::uint32_t rec_idx_a, std::uint32_t rec_idx_b, float dist, const ContactContext& ctx) -> InteractionType {
//       const auto& params = ctx.get_params<HalogenParams>();
//       const auto &donor    = ctx.topology.atom(rec_idx_a).atom;
//       const auto &acceptor = ctx.topology.atom(rec_idx_b).atom;
// 
//       if (!halo_geo::are_geometrically_viable(ctx.molecule(), donor, acceptor, params)) return InteractionType::None;
// 
//       return InteractionType::Halogen;
//     }
//   );
// }

} // namespace lahuta
