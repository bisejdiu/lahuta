#include "contacts/halogen_bonds.hpp"
#include "contacts/geometry.hpp"
#include "contacts/search.hpp"
#include "lahuta.hpp"

namespace lahuta {

AtomType add_halogen_donor(const RDKit::RWMol &mol, const RDKit::Atom &atom) {
  if (HalogenDonors.count(atom.getAtomicNum())) return AtomType::XBOND_DONOR;
  return AtomType::NONE;
}

AtomType add_halogen_acceptor(const RDKit::RWMol &mol, const RDKit::Atom &atom) {
  if (!HalogenAcceptors.count(atom.getAtomicNum())) return AtomType::NONE;

  for (const auto &bond : mol.atomBonds(&atom)) {
    const RDKit::Atom *neighbor = bond->getOtherAtom(&atom);
    int neighbor_atomic_num = neighbor->getAtomicNum();

    if (HalogenBinders.count(neighbor_atomic_num)) return AtomType::XBOND_ACCEPTOR;
  }

  return AtomType::NONE;
}

bool are_geometrically_viable(
    const RDKit::RWMol &mol, const RDKit::Atom &donor, const RDKit::Atom &acceptor,
    const HalogenParams &opts) {

  auto [halogen_angles, _] = geometry::calculate_angle(mol, donor, acceptor, true);
  if (halogen_angles.size() != 1) return false;

  if (opts.optimal_angle - halogen_angles[0] > opts.angle_max) return false;

  auto [acceptor_angles, __] = geometry::calculate_angle(mol, acceptor, donor, true);
  if (acceptor_angles.empty()) return false;

  bool exit_outer_flag = false;
  for (double acceptor_angle : acceptor_angles) {
    if (opts.optimal_acceptor_angle - acceptor_angle > opts.angle_max) {
      exit_outer_flag = true;
      continue;
    }
  }
  if (exit_outer_flag) return false;

  return true;
}

Contacts find_halogen_bonds(const Luni &luni, HalogenParams opts) {

  Contacts contacts(&luni);

  const auto donor_atoms = AtomEntityCollection::filter(&luni, AtomType::XBOND_DONOR);
  const auto acceptor_atoms = AtomEntityCollection::filter(&luni, AtomType::XBOND_ACCEPTOR);

  auto nbrs = EntityNeighborSearch::search(donor_atoms, acceptor_atoms, opts.distance_max);

  for (const auto &[pair, dist] : nbrs) {
    auto [donor_index, acceptor_index] = pair;
    const auto &donor = donor_atoms.get_data()[donor_index];
    const auto &acceptor = acceptor_atoms.get_data()[acceptor_index];

    if (!are_geometrically_viable(luni.get_molecule(), *donor.atom, *acceptor.atom, opts)) continue;

    contacts.add(Contact(
        static_cast<EntityID>(donor.atom->getIdx()),
        static_cast<EntityID>(acceptor.atom->getIdx()),
        dist,
        InteractionType::Halogen));
  }

  return contacts;
}

} // namespace lahuta
