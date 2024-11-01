#include "contacts/hydrogen_bonds.hpp"
#include "contacts/geometry.hpp"
#include "contacts/search.hpp"
#include "contacts/utils.hpp"
#include "lahuta.hpp"

namespace lahuta {

bool is_water(const RDKit::Atom &atom) {
  auto res_info = static_cast<const RDKit::AtomPDBResidueInfo *>(atom.getMonomerInfo());

  if (res_info) {
    return std::find(WaterResNames.begin(), WaterResNames.end(), res_info->getResidueName())
           != WaterResNames.end();
  }
  return false;
}

auto *closest_hydrogen_atom(const RDKit::RWMol &mol, const RDKit::Atom &atom_a, const RDKit::Atom &atom_b) {
  auto &conf = mol.getConformer();
  auto &pos_a = conf.getAtomPos(atom_a.getIdx());
  auto &pos_b = conf.getAtomPos(atom_b.getIdx());

  double min_dist_sq = (pos_a - pos_b).lengthSq();
  auto *closest_hydrogen = &atom_a; // Initialize to atom_a itself

  for (const auto bond : mol.atomBonds(&atom_a)) {
    int other_atom_idx = bond->getOtherAtomIdx(atom_a.getIdx());
    auto *neighbor = mol.getAtomWithIdx(other_atom_idx);
    if (neighbor->getAtomicNum() == 1) { // Hydrogen
      auto &pos_h = conf.getAtomPos(neighbor->getIdx());
      double dist_sq = (pos_h - pos_b).lengthSq();
      if (dist_sq < min_dist_sq) {
        min_dist_sq = dist_sq;
        closest_hydrogen = neighbor;
      }
    }
  }
  return closest_hydrogen;
}

std::vector<const RDGeom::Point3D *> get_neighbor_positions(
    const RDKit::Atom &atom_a, const RDKit::Conformer &conf, const RDKit::RWMol &mol, bool ignore_hydrogens) {
  std::vector<const RDGeom::Point3D *> neighbor_positions;
  const RDKit::Atom *first_neighbor = nullptr;

  // First pass: try to find direct neighbors
  for (const auto &bond : mol.atomBonds(&atom_a)) {
    auto *neighbor = bond->getOtherAtom(&atom_a);
    if (ignore_hydrogens && neighbor->getAtomicNum() == 1) {
      continue;
    }

    const RDGeom::Point3D &neighbor_pos = conf.getAtomPos(neighbor->getIdx());
    neighbor_positions.push_back(&neighbor_pos);

    if (neighbor_positions.size() == 1) {
      first_neighbor = neighbor;
    } else if (neighbor_positions.size() == 2) {
      break;
    }
  }

  // If only one neighbor was found, look at its neighbors
  if (neighbor_positions.size() == 1 && first_neighbor != nullptr) {
    for (const auto &bond : mol.atomBonds(first_neighbor)) {
      auto *next_neighbor = bond->getOtherAtom(first_neighbor);
      if (next_neighbor == &atom_a || (ignore_hydrogens && next_neighbor->getAtomicNum() == 1)) {
        continue;
      }

      const RDGeom::Point3D &neighbor_pos = conf.getAtomPos(next_neighbor->getIdx());
      neighbor_positions.push_back(&neighbor_pos);
      break;
    }
  }

  return neighbor_positions;
}

bool in_aromatic_ring_with_N_or_O(const RDKit::RWMol &mol, const RDKit::Atom &atom) {
  unsigned int atomIdx = atom.getIdx();
  const RDKit::RingInfo *ringInfo = mol.getRingInfo();

  for (const auto &ring : ringInfo->atomRings()) {
    if (std::find(ring.begin(), ring.end(), atomIdx) != ring.end()) {
      // Check if the ring is aromatic
      bool isAromatic = true;
      for (unsigned int ringAtomIdx : ring) {
        if (!mol.getAtomWithIdx(ringAtomIdx)->getIsAromatic()) {
          isAromatic = false;
          break;
        }
      }
      if (!isAromatic) {
        continue;
      }

      // Check if the ring contains an electronegative element (N or O)
      for (unsigned int ringAtomIdx : ring) {
        const RDKit::Atom *ringAtom = mol.getAtomWithIdx(ringAtomIdx);
        int atomicNum = ringAtom->getAtomicNum();
        if (atomicNum == 7 || atomicNum == 8) { // N or O
          // electronegative element (N or O) in the aromatic ring
          return true;
        }
      }
    }
  }

  return false;
}

bool are_geometrically_viable(
    const RDKit::RWMol &mol, const RDKit::Atom &donor, const RDKit::Atom &acceptor,
    const GeometryOptions &opts) {

  // donor angles
  const auto &[don_angles, don_h_angles] = calculate_angle(mol, donor, acceptor, opts.ignore_hydrogens);

  double ideal_don_angle = get_atom_geometry_angle(donor.getHybridization());

  for (double don_angle : don_angles) {
    if (std::abs(ideal_don_angle - don_angle) > opts.max_don_angle_dev) {
      return false;
    }
  }

  if (!don_h_angles.empty() && std::all_of(don_h_angles.begin(), don_h_angles.end(), [&opts](double h_angle) {
        return h_angle >= opts.max_don_angle_dev;
      })) {
    return false;
  }

  // out-of-plane angle for sp2 hybridized atoms
  if (donor.getHybridization() == HybridizationType::SP2) {
    double out_of_plane = compute_plane_angle(mol, donor, acceptor);
    if (out_of_plane > opts.max_don_out_of_plane_angle) {
      return false;
    }
  }

  // Potentially update donor index to the closest hydrogen
  const RDKit::Atom *donor_atom_for_acc = &donor;
  if (!opts.ignore_hydrogens && !don_h_angles.empty()) {
    donor_atom_for_acc = closest_hydrogen_atom(mol, donor, acceptor);
  }

  // acceptor angles
  const auto &[acc_angles, acc_h_angles] =
      calculate_angle(mol, acceptor, *donor_atom_for_acc, opts.ignore_hydrogens);

  double ideal_acc_angle = get_atom_geometry_angle(acceptor.getHybridization());

  // Check acceptor angles (do not limit large acceptor angles)
  for (double acc_angle : acc_angles) {
    if (ideal_acc_angle - acc_angle > opts.max_acc_angle_dev) {
      return false;
    }
  }
  for (double acc_h_angle : acc_h_angles) {
    if (ideal_acc_angle - acc_h_angle > opts.max_acc_angle_dev) {
      return false;
    }
  }

  // out-of-plane angle for sp2 hybridized atoms
  if (acceptor.getHybridization() == RDKit::Atom::HybridizationType::SP2) {
    double out_of_plane = compute_plane_angle(mol, acceptor, *donor_atom_for_acc);
    if (out_of_plane > opts.max_acc_out_of_plane_angle) {
      return false;
    }
  }
  return true;
}

AtomType add_hydrogen_donor(const RDKit::RWMol &mol, const RDKit::Atom &atom) {
  int implicit_h = atom.getProp<int>("computed_implicit_h");
  int total_h = atom.getNumExplicitHs() + implicit_h;

  // include both nitrogen atoms in histidine due to
  // their often ambiguous protonation assignment
  if (is_histidine_nitrogen(atom, mol)) {
    return AtomType::HBOND_DONOR;
  }

  // nitrogen, oxygen, or sulfur with hydrogen attached
  int atomic_num = atom.getAtomicNum();
  if (total_h > 0 && (atomic_num == 7 || atomic_num == 8 || atomic_num == 16)) {
    return AtomType::HBOND_DONOR;
  }

  return AtomType::NONE;
}

AtomType add_hydrogen_acceptor(const RDKit::RWMol &mol, const RDKit::Atom &atom) {
  auto _info = atom.getMonomerInfo();
  auto *res_info = static_cast<const RDKit::AtomPDBResidueInfo *>(_info);

  int implicit_h = atom.getProp<int>("computed_implicit_h");
  int formal_charge = atom.getFormalCharge();
  auto hybridization = atom.getHybridization();
  int atomic_num = atom.getAtomicNum();

  // Assume all oxygen atoms are acceptors!
  if (atomic_num == 8) {
    return AtomType::HBOND_ACCEPTOR;
  }

  else if (atomic_num == 7) {
    if (is_histidine_nitrogen(atom, mol)) {
      // include both nitrogen atoms in histidine due to
      // their often ambiguous protonation assignment
      return AtomType::HBOND_ACCEPTOR;
    } else if (formal_charge < 1) {
      // Neutral nitrogen might be an acceptor
      // It must have at least one lone pair not conjugated
      unsigned int total_bonds = get_bond_count(mol, atom) + implicit_h;

      if ((hybridization == HybridizationType::SP3 && total_bonds < 4)
          || (hybridization == HybridizationType::SP2 && total_bonds < 3)
          || (hybridization == HybridizationType::SP && total_bonds < 2)) {
        return AtomType::HBOND_ACCEPTOR;
      }
    }
  }

  else if (atomic_num == 16) {
    if (res_info) {
      // FIX: hardcoded residue names
      std::string res_name = res_info->getResidueName();
      if (res_name == "CYS" || res_name == "MET" || formal_charge == -1) {
        return AtomType::HBOND_ACCEPTOR;
      }
    }
  }

  return AtomType::NONE;
}

AtomType add_weak_hydrogen_donor(const RDKit::RWMol &mol, const RDKit::Atom &atom) {
  int implicit_h = atom.getProp<int>("computed_implicit_h");
  int total_h = atom.getNumExplicitHs() + implicit_h;
  auto info = static_cast<const RDKit::AtomPDBResidueInfo *>(atom.getMonomerInfo());
  /*std::cout << "Debug: " << atom.getIdx() << " " */
  /*  << info->getResidueName() << " " << info->getName() << " "*/
  /*  << implicit_h << " " << total_h << std::endl;*/

  // Check if the atom is carbon
  if (atom.getAtomicNum() == 6 && total_h > 0) {
    // Check if the atom is bonded to nitrogen or oxygen, or in an aromatic ring
    // with N or O
    if (get_bond_count(mol, atom, 7) > 0 ||        // Bonded to nitrogen
        get_bond_count(mol, atom, 8) > 0 ||        // Bonded to oxygen
        in_aromatic_ring_with_N_or_O(mol, atom)) { // Aromatic ring with N or O
      return AtomType::WEAK_HBOND_DONOR;
    }
  }

  return AtomType::NONE;
}

void find_hydrogen_bonds(Luni &luni, const GeometryOptions &opts, Contacts &container) {

  const auto &mol = luni.get_molecule();

  const AtomDataVec donor_atoms = get_atom_data(&luni, AtomType::HBOND_DONOR);
  const AtomDataVec acceptor_atoms = get_atom_data(&luni, AtomType::HBOND_ACCEPTOR);

  auto max_dist = std::max(opts.max_dist_sq, opts.max_sulfur_dist_sq);

  EntityNeighborSearch ens(mol.getConformer());
  auto results = ens.search(donor_atoms, acceptor_atoms, std::sqrt(max_dist));

  for (const auto &[pair, dist] : results) {
    auto [donor_index, acceptor_index] = pair;
    const auto &donor = donor_atoms.get_data()[donor_index];
    const auto &acceptor = acceptor_atoms.get_data()[acceptor_index];

    // FIX: improve this (it only checks residue numbers)
    if (are_residueids_close(mol, *donor.atom, *acceptor.atom, 1)) {
      continue;
    }

    // TODO: Potential optimization by searching only around S atoms.
    // Would likely complicate the code with only minor performance gains
    double max_dist_sq;
    if (donor.atom->getAtomicNum() == 16 || acceptor.atom->getAtomicNum() == 16) {
      max_dist_sq = opts.max_sulfur_dist_sq;
    } else {
      max_dist_sq = opts.max_dist_sq;
    }
    if (dist > max_dist_sq) {
      continue;
    }

    if (!opts.include_water && is_water_hbond(*donor.atom, *acceptor.atom)) {
      continue;
    }

    if (!are_geometrically_viable(mol, *donor.atom, *acceptor.atom, opts)) {
      continue;
    }

    container.add(Contact(
        static_cast<EntityID>(donor.atom->getIdx()),
        static_cast<EntityID>(acceptor.atom->getIdx()),
        dist,
        InteractionType::HydrogenBond));
  }
}

void find_weak_hydrogen_bonds(Luni &luni, const GeometryOptions &opts, Contacts &container) {

  const auto &mol = luni.get_molecule();

  const auto weak_donor_atoms = get_atom_data(&luni, AtomType::WEAK_HBOND_DONOR);
  const auto acceptor_atoms = get_atom_data(&luni, AtomType::HBOND_ACCEPTOR);

  EntityNeighborSearch ens(mol.getConformer());
  auto nbrs = ens.search(weak_donor_atoms, acceptor_atoms, std::sqrt(opts.max_dist_sq));

  for (const auto &[pair, dist] : nbrs) {
    auto [wdonor_index, acceptor_index] = pair;
    const auto &wdonor = weak_donor_atoms.get_data()[wdonor_index];
    const auto &acceptor = acceptor_atoms.get_data()[acceptor_index];

    if (are_residueids_close(mol, *wdonor.atom, *acceptor.atom, 1)) {
      continue;
    }

    if (!opts.include_water && is_water_hbond(*wdonor.atom, *acceptor.atom)) {
      continue;
    }

    if (!are_geometrically_viable(mol, *wdonor.atom, *acceptor.atom, opts)) {
      continue;
    }

    container.add(Contact(
        static_cast<EntityID>(wdonor.atom->getIdx()),
        static_cast<EntityID>(acceptor.atom->getIdx()),
        dist,
        InteractionType::WeakHydrogenBond));
  }
}

} // namespace lahuta
