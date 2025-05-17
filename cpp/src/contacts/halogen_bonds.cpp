#include "contacts/halogen_bonds.hpp"
#include "contacts/halo_geo_validity.hpp"
#include "contacts/search.hpp"
#include "lahuta.hpp"
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

Contacts find_halogen_bonds(const Luni &luni, std::optional<HalogenParams> params) {

  Contacts contacts(&luni);
  HalogenParams opts = params.value_or(HalogenParams{});

  const auto donor_atoms    = AtomEntityCollection::filter(&luni, AtomType::XBOND_DONOR);
  const auto acceptor_atoms = AtomEntityCollection::filter(&luni, AtomType::XBOND_ACCEPTOR);

  auto nbrs = EntityNeighborSearch::search(donor_atoms, acceptor_atoms, opts.distance_max);

  for (const auto &[pair, dist] : nbrs) {
    auto [donor_index, acceptor_index] = pair;
    const auto &donor    = donor_atoms   .get_data()[donor_index];
    const auto &acceptor = acceptor_atoms.get_data()[acceptor_index];

    if (!halo_geo::are_geometrically_viable(luni.get_molecule(), *donor.atom, *acceptor.atom, opts)) continue;

    contacts.add(Contact(
        static_cast<EntityID>(donor.atom->getIdx()),
        static_cast<EntityID>(acceptor.atom->getIdx()),
        dist,
        InteractionType::Halogen));
  }

  return contacts;
}

} // namespace lahuta
