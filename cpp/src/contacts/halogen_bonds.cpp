#include "contacts/halogen_bonds.hpp"
#include "contacts/geometry.hpp"
#include "contacts/search.hpp"
#include "lahuta.hpp"

namespace lahuta {

AtomType add_halogen_donor(const RDKit::RWMol &mol, const RDKit::Atom &atom) {
  if (hal_bond_elements.count(atom.getAtomicNum())) {
    return AtomType::XBOND_DONOR;
  }
  return AtomType::NONE;
}

AtomType add_halogen_acceptor(const RDKit::RWMol &mol, const RDKit::Atom &atom) {
  if (X.count(atom.getAtomicNum())) {
    bool flag = false;

    for (const auto &bond : mol.atomBonds(&atom)) {
      const RDKit::Atom *neighbor = bond->getOtherAtom(&atom);
      int neighbor_atomic_num = neighbor->getAtomicNum();

      if (Y.count(neighbor_atomic_num)) {
        flag = true;
        break;
      }
    }

    if (flag) {
      return AtomType::XBOND_ACCEPTOR;
    }
  }
  return AtomType::NONE;
}

Contacts find_halogen_bonds(const Luni &luni, HalogenParams opts) {
  // NOTE: we do not check for same residue here

  Contacts contacts(&luni);
  const auto &mol = luni.get_molecule();

  const auto donor_atoms = get_atom_data(&luni, AtomType::XBOND_DONOR);
  const auto acceptor_atoms = get_atom_data(&luni, AtomType::XBOND_ACCEPTOR);

  EntityNeighborSearch ens(mol.getConformer());
  auto nbrs = ens.search(donor_atoms, acceptor_atoms, opts.distance_max);

  for (const auto &[pair, dist] : nbrs) {

    auto [donor_index, acceptor_index] = pair;
    const auto &donor = donor_atoms.get_data()[donor_index];
    const auto &acceptor = acceptor_atoms.get_data()[acceptor_index];

    // TODO: put this into a `validate_geometry` function
    auto [halogen_angles, _] = calculate_angle(mol, *donor.atom, *acceptor.atom, true);
    if (halogen_angles.size() != 1) {
      continue;
    }
    if (opts.optimal_angle - halogen_angles[0] > opts.angle_max) continue;

    auto [acceptor_angles, __] = calculate_angle(mol, *acceptor.atom, *donor.atom, true);
    if (acceptor_angles.empty()) continue;

    bool exit_outer_flag = false;
    for (double acceptor_angle : acceptor_angles) {
      if (opts.optimal_acceptor_angle - acceptor_angle > opts.angle_max) {
        exit_outer_flag = true;
        continue;
      }
    }
    if (exit_outer_flag) continue;

    contacts.add(Contact(
        static_cast<EntityID>(donor.atom->getIdx()),
        static_cast<EntityID>(acceptor.atom->getIdx()),
        dist,
        InteractionType::Halogen));
  }

  return contacts;
}

} // namespace lahuta
