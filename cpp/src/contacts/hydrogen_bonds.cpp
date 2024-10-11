#include "contacts/hydrogen_bonds.hpp"
#include "contacts/geometry.hpp"

namespace lahuta {

bool is_water(const RDKit::Atom &atom) {
  auto res_info =
      static_cast<const RDKit::AtomPDBResidueInfo *>(atom.getMonomerInfo());

  if (res_info) {
    return std::find(WaterResNames.begin(), WaterResNames.end(),
                     res_info->getResidueName()) != WaterResNames.end();
  }
  return false;
}

auto *closest_hydrogen_atom(const RDKit::RWMol &mol, const RDKit::Atom &atom_a,
                            const RDKit::Atom &atom_b) {
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

std::vector<const RDGeom::Point3D *>
get_neighbor_positions(const RDKit::Atom &atom_a, const RDKit::Conformer &conf,
                       const RDKit::RWMol &mol, bool ignore_hydrogens) {
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
      if (next_neighbor == &atom_a ||
          (ignore_hydrogens && next_neighbor->getAtomicNum() == 1)) {
        continue;
      }

      const RDGeom::Point3D &neighbor_pos =
          conf.getAtomPos(next_neighbor->getIdx());
      neighbor_positions.push_back(&neighbor_pos);
      break;
    }
  }

  return neighbor_positions;
}

bool check_geometry_constraints(const RDKit::RWMol &mol,
                                const RDKit::Atom &donor,
                                const RDKit::Atom &acceptor,
                                const GeometryOptions &opts) {

  // donor angles
  const auto &[don_angles, don_h_angles] =
      calculate_angle(mol, donor, acceptor, opts.ignore_hydrogens);

  double ideal_don_angle = get_atom_geometry_angle(donor.getHybridization());

  for (double don_angle : don_angles) {
    if (std::abs(ideal_don_angle - don_angle) > opts.max_don_angle_dev) {
      return false;
    }
  }

  if (!don_h_angles.empty() &&
      std::all_of(don_h_angles.begin(), don_h_angles.end(),
                  [&opts](double h_angle) {
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
  const auto &[acc_angles, acc_h_angles] = calculate_angle(
      mol, acceptor, *donor_atom_for_acc, opts.ignore_hydrogens);

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
    double out_of_plane =
        compute_plane_angle(mol, acceptor, *donor_atom_for_acc);
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

AtomType add_hydrogen_acceptor(const RDKit::RWMol &mol,
                               const RDKit::Atom &atom) {
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

      if ((hybridization == HybridizationType::SP3 && total_bonds < 4) ||
          (hybridization == HybridizationType::SP2 && total_bonds < 3) ||
          (hybridization == HybridizationType::SP && total_bonds < 2)) {
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

} // namespace lahuta
