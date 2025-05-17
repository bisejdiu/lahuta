#include "contacts/halogen_bonds.hpp"
#include "contacts/halo_geo_validity.hpp"
#include "contacts/search.hpp"
#include "lahuta.hpp"

namespace lahuta {

std::unordered_set<int> HalogenDonors    = {17, 35, 53};
std::unordered_set<int> HalogenAcceptors = {7, 8, 16};
std::unordered_set<int> HalogenBinders   = {6, 7, 15, 16};

AtomType add_halogen_donor(const RDKit::RWMol &mol, const RDKit::Atom &atom) {
  if (HalogenDonors.count(atom.getAtomicNum())) return AtomType::XBOND_DONOR;
  return AtomType::NONE;
}

AtomType add_halogen_acceptor(const RDKit::RWMol &mol, const RDKit::Atom &atom) {
  if (!HalogenAcceptors.count(atom.getAtomicNum())) return AtomType::NONE;

  for (const auto bond : mol.atomBonds(&atom)) {
    const RDKit::Atom *neighbor = bond->getOtherAtom(&atom);
    int neighbor_atomic_num = neighbor->getAtomicNum();

    if (HalogenBinders.count(neighbor_atomic_num)) return AtomType::XBOND_ACCEPTOR;
  }

  return AtomType::NONE;
}

Contacts find_halogen_bonds(const Luni &luni, HalogenParams opts) {

  Contacts contacts(&luni);

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
