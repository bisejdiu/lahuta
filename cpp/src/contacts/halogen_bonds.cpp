#include "contacts/halogen_bonds.hpp"
#include "contacts/halo_geo_validity.hpp"
#include "entities/find_contacts.hpp"
#include <elements.hpp>

namespace lahuta {

std::unordered_set<Element> HalogenDonors    = {Element::Cl, Element::Br, Element::I};
std::unordered_set<Element> HalogenAcceptors = {Element::N,  Element::O,  Element::S};
std::unordered_set<Element> HalogenBinders   = {Element::C,  Element::N,  Element::P, Element::S};

AtomType add_halogen_donor(const RDKit::RWMol &mol, const RDKit::Atom &atom) {
  const auto at_n = static_cast<Element>(atom.getAtomicNum());
  if (HalogenDonors.count(at_n)) return AtomType::XBOND_DONOR;
  return AtomType::NONE;
}

AtomType add_halogen_acceptor(const RDKit::RWMol &mol, const RDKit::Atom &atom) {
  const auto at_n = static_cast<Element>(atom.getAtomicNum());
  if (!HalogenAcceptors.count(at_n)) return AtomType::NONE;

  for (const auto bond : mol.atomBonds(&atom)) {
    const RDKit::Atom *nbr = bond->getOtherAtom(&atom);
    const auto nbr_at_n = static_cast<Element>(nbr->getAtomicNum());

    if (HalogenBinders.count(nbr_at_n)) return AtomType::XBOND_ACCEPTOR;
  }

  return AtomType::NONE;
}

ContactSet find_halogen_bonds(const Topology &topology, const HalogenParams &params) {
  return find_contacts(
    topology,
    [](const AtomRec &rec) { return (rec.type & AtomType::XBOND_DONOR)    == AtomType::XBOND_DONOR; },
    [](const AtomRec &rec) { return (rec.type & AtomType::XBOND_ACCEPTOR) == AtomType::XBOND_ACCEPTOR; },
    {params.distance_max, 0, 0, 0.7},
    [&topology, &params](std::uint32_t rec_idx_a, std::uint32_t rec_idx_b, float dist) -> InteractionType {
      const auto &donor_rec    = topology.atom(rec_idx_a);
      const auto &acceptor_rec = topology.atom(rec_idx_b);

      const auto &mol = topology.molecule();
      const auto *donor_atom    = mol.getAtomWithIdx(donor_rec.idx);
      const auto *acceptor_atom = mol.getAtomWithIdx(acceptor_rec.idx);

      if (!halo_geo::are_geometrically_viable(mol, *donor_atom, *acceptor_atom, params)) return InteractionType::None;

      return InteractionType::Halogen;
    }
  );
}

} // namespace lahuta
