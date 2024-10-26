#include "contacts/halogen_bonds.hpp"
#include "contacts/geometry.hpp"
#include "contacts/utils.hpp"
#include "lahuta.hpp"

namespace lahuta {

const double OptimalHalogenAngle = M_PI;           // 180 degrees in radians
const double OptimalAcceptorAngle = 2.09439510239; // 120 degrees in radians

// Convert options (HalogenBondsProps) to HalogenBondsOptions
HalogenBondsOptions get_options(const HalogenBondsParams &params) {
  return {deg_to_rad(params.angleMax)}; //
}

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

void find_halogen_bonds(Luni &luni, const GeometryOptions &opts, Contacts &container) {
  // NOTE: we do not check for same residue here

  const auto &mol = luni.get_molecule();

  const auto donor_atoms = get_atom_data(&luni, AtomType::XBOND_DONOR);
  const auto acceptor_atoms = get_atom_data(&luni, AtomType::XBOND_ACCEPTOR);

  double max_dist_sq = 4.0 * 4.0;

  auto grid = FastNS(acceptor_atoms.positions(), std::sqrt(max_dist_sq));
  auto nbrs = grid.search(donor_atoms.positions());

  // Halogen bond options (currently using default parameters)
  HalogenBondsOptions opts_h = get_options(HalogenBondsParams{});

  for (const auto &[pair, dist] : nbrs) {

    auto [donor_index, acceptor_index] = pair;
    const auto &donor = donor_atoms.data[donor_index];
    const auto &acceptor = acceptor_atoms.data[acceptor_index];

    // TODO: put this into a `validate_geometry` function
    auto [halogen_angles, _] = calculate_angle(mol, *donor.atom, *acceptor.atom, true);
    if (halogen_angles.size() != 1) {
      continue;
    }
    if (OptimalHalogenAngle - halogen_angles[0] > opts_h.angleMax) {
      continue;
    }

    auto [acceptor_angles, __] = calculate_angle(mol, *acceptor.atom, *donor.atom, true);
    if (acceptor_angles.empty()) {
      continue;
    }
    bool exit_outer_flag = false;
    for (double acceptor_angle : acceptor_angles) {
      if (OptimalAcceptorAngle - acceptor_angle > opts_h.angleMax) {
        exit_outer_flag = true;
        continue;
      }
    }
    if (exit_outer_flag) {
      continue;
    }

    container.add(Contact(
        EntityID(donor.atom->getIdx()), EntityID(acceptor.atom->getIdx()), dist, InteractionType::Halogen));
  }
}

} // namespace lahuta
