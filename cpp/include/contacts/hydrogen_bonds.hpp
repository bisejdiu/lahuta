/*
    This code adapts functionality from the Mol* (mol-star) molecular visualization package,
    originally in TypeScript, for use in C++ with RDKit in the Lahuta project.

    Mol* is licensed under the MIT License (https://github.com/molstar/molstar),
    while this adaptation is released under the GNU General Public License (GPL).
*/
#ifndef LAHUTA_HYDROGEN_BONDS_HPP
#define LAHUTA_HYDROGEN_BONDS_HPP

#include "atom_types.hpp"
#include "nn.hpp"
#include "utils.hpp"
#include "valence_model.hpp"
#include <GraphMol/Atom.h>
#include <GraphMol/MonomerInfo.h>
#include <GraphMol/RWMol.h>

using HybridizationType = RDKit::Atom::HybridizationType;

namespace lahuta {

class Luni;

constexpr double MAX_DIST = 3.5;
constexpr double MAX_SULFUR_DIST = 4.1;
constexpr double MAX_ACC_ANGLE_DEV = M_PI / 4.0;          // 45 degrees
constexpr double MAX_DON_ANGLE_DEV = M_PI / 4.0;          // 45 degrees
constexpr double MAX_DON_OUT_OF_PLANE_ANGLE = M_PI / 4.0; // 45 degrees
constexpr double MAX_ACC_OUT_OF_PLANE_ANGLE = M_PI / 2.0; // 90 degrees

inline struct GeometryOptions {
  /// Ignore hydrogens in geometry calculations
  bool ignore_hydrogens = false;
  /// Include backbone-to-backbone hydrogen bonds
  bool include_backbone = true;
  /// Include water-to-water hydrogen bonds
  bool include_water = true;

  /// Maximum distance for sulfur atoms
  double max_sulfur_dist = MAX_SULFUR_DIST;
  /// Maximum distance for hydrogen bonds
  double max_dist = MAX_DIST;

  /// Maximum deviation from ideal acceptor angle
  double max_acc_angle_dev = MAX_ACC_ANGLE_DEV;
  /// Maximum deviation from ideal donor angle
  double max_don_angle_dev = MAX_DON_ANGLE_DEV;
  /// Maximum out-of-plane deviation for acceptor
  double max_acc_out_of_plane_angle = MAX_ACC_OUT_OF_PLANE_ANGLE;
  /// Maximum out-of-plane deviation for donor
  double max_don_out_of_plane_angle = MAX_DON_OUT_OF_PLANE_ANGLE;

} hydrogen_bond_opts;

constexpr double get_atom_geometry_angle(HybridizationType hybridization) {
  switch (hybridization) {
    case HybridizationType::SP: // 180 degrees
      return M_PI;
    case HybridizationType::SP2: // 120 degrees
      return M_PI * 2.0 / 3.0;
    case HybridizationType::SP3:
      return deg_to_rad(109.4721);
    case HybridizationType::SP3D2: // 90 degrees
      return M_PI / 2.0;

    default:
      return M_PI * 2.0 / 3.0;
  }
}

// Residue names that represent water molecules
const std::array<std::string, 11> WaterResNames = {
    "SOL", "WAT", "HOH", "H2O", "W", "DOD", "D3O", "TIP", "TIP3", "TIP4", "SPC"};

const std::vector<std::string> HistidineResNames = {"HIS", "HID", "HIE", "HIP"};

bool is_water(const RDKit::Atom &atom);

inline bool is_water_hbond(const RDKit::Atom &atom_a, const RDKit::Atom &atom_b) {
  return is_water(atom_a) && is_water(atom_b);
}

// FIX: Add support for other histidine resname variants
inline bool is_histidine_nitrogen(const RDKit::Atom &atom, const RDKit::RWMol &mol) {
  auto res_info = static_cast<const RDKit::AtomPDBResidueInfo *>(atom.getMonomerInfo());

  // FIX: This requires the initialization of ringInfo
  // NOTE: should be fixed now. 
  auto ri = mol.getRingInfo();
  bool is_in_ring = ri->numAtomRings(atom.getIdx()) > 0;
  if (res_info && res_info->getResidueName() == "HIS" && atom.getAtomicNum() == 7 && is_in_ring) {
    return true;
  }
  return false;
}

std::vector<const RDGeom::Point3D *> get_neighbor_positions(
    const RDKit::Atom &atom_a, const RDKit::Conformer &conf, const RDKit::RWMol &mol,
    bool ignore_hydrogens = true);

// Find the closest hydrogen atom bonded to atom_a with respect to atom_b
auto *closest_hydrogen_atom(const RDKit::RWMol &mol, const RDKit::Atom &atom_a, const RDKit::Atom &atom_b);

// Validate geometry for hydrogen bonds
bool are_geometrically_viable(
    const RDKit::RWMol &mol, const RDKit::Atom &donor, const RDKit::Atom &acceptor,
    const GeometryOptions &opts);

/// atom is in an aromatic ring with an electronegative element (N or O)
bool in_aromatic_ring_with_N_or_O(const RDKit::RWMol &mol, const RDKit::Atom &atom);

AtomType add_hydrogen_donor(const RDKit::RWMol &mol, const RDKit::Atom &atom);
AtomType add_hydrogen_acceptor(const RDKit::RWMol &mol, const RDKit::Atom &atom);
AtomType add_weak_hydrogen_donor(const RDKit::RWMol &mol, const RDKit::Atom &atom);

Contacts find_hydrogen_bonds(const Luni &luni, const GeometryOptions &opts = hydrogen_bond_opts);
Contacts find_weak_hydrogen_bonds(const Luni &luni, const GeometryOptions &opts = hydrogen_bond_opts);

} // namespace lahuta

#endif // LAHUTA_HYDROGEN_BONDS_HPP
