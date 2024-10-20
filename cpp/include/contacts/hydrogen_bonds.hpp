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
#include "nsgrid.hpp"
#include "valence_model.hpp"
#include <GraphMol/Atom.h>
#include <GraphMol/MonomerInfo.h>
#include <GraphMol/RWMol.h>

using HybridizationType = RDKit::Atom::HybridizationType;

constexpr float MAX_DIST_SQ = 12.25;
constexpr float MAX_SULFUR_DIST_SQ = 16.81;
constexpr float MAX_ACC_ANGLE_DEV = 45.0;
constexpr float MAX_DON_ANGLE_DEV = 45.0;
constexpr float MAX_DON_OUT_OF_PLANE_ANGLE = 45.0;
constexpr float MAX_ACC_OUT_OF_PLANE_ANGLE = 90.0;

namespace lahuta {

class Luni;

// Residue names that represent water molecules
const std::array<std::string, 11> WaterResNames = {
    "SOL", "WAT", "HOH", "H2O", "W", "DOD", "D3O", "TIP", "TIP3", "TIP4", "SPC"};

const std::vector<std::string> HistidineResNames = {"HIS", "HID", "HIE", "HIP"};

inline double deg_to_rad(double degrees) { return degrees * (M_PI / 180.0); }

bool is_water(const RDKit::Atom &atom);

inline bool is_water_hbond(const RDKit::Atom &atom_a, const RDKit::Atom &atom_b) {
  return is_water(atom_a) && is_water(atom_b);
}

// FIX: Add support for other histidine resname variants
inline bool is_histidine_nitrogen(const RDKit::Atom &atom, const RDKit::RWMol &mol) {
  auto res_info = static_cast<const RDKit::AtomPDBResidueInfo *>(atom.getMonomerInfo());

  // FIX: This requires the initialization of ringInfo
  auto ri = mol.getRingInfo();
  bool is_in_ring = ri->numAtomRings(atom.getIdx()) > 0;
  if (res_info && res_info->getResidueName() == "HIS" && atom.getAtomicNum() == 7 && is_in_ring) {
    return true;
  }
  return false;
}

struct GeometryOptions {
  // Ignore hydrogens in geometry calculations
  bool ignore_hydrogens;
  // Include backbone-to-backbone hydrogen bonds
  bool include_backbone;
  // Maximum deviation from ideal acceptor angle
  double max_acc_angle_dev;
  // Maximum deviation from ideal donor angle
  double max_don_angle_dev;
  // Maximum out-of-plane deviation for acceptor
  double max_acc_out_of_plane_angle;
  // Maximum out-of-plane deviation for donor
  double max_don_out_of_plane_angle;
  // Include water-to-water hydrogen bonds
  bool include_water;
  // Maximum distance for sulfur atoms
  double max_sulfur_dist_sq;
  // Maximum distance for hydrogen bonds
  double max_dist_sq;

  GeometryOptions(
      bool _ignore_hydrogens = false, bool _include_backbone = true,
      double _max_acc_angle_dev = deg_to_rad(MAX_ACC_ANGLE_DEV),
      double _max_don_angle_dev = deg_to_rad(MAX_DON_ANGLE_DEV),
      double _max_acc_out_of_plane_angle = deg_to_rad(MAX_ACC_OUT_OF_PLANE_ANGLE),
      double _max_don_out_of_plane_angle = deg_to_rad(MAX_DON_OUT_OF_PLANE_ANGLE), bool _include_water = true,
      double _max_sulfur_dist_sq = MAX_SULFUR_DIST_SQ, double _max_dist_sq = MAX_DIST_SQ)
      : ignore_hydrogens(_ignore_hydrogens), include_backbone(_include_backbone),
        max_acc_angle_dev(_max_acc_angle_dev), max_don_angle_dev(_max_don_angle_dev),
        max_acc_out_of_plane_angle(_max_acc_out_of_plane_angle),
        max_don_out_of_plane_angle(_max_don_out_of_plane_angle), include_water(_include_water),
        max_sulfur_dist_sq(_max_sulfur_dist_sq), max_dist_sq(_max_dist_sq) {}

  void set_max_dist(double max_dist) { max_dist_sq = max_dist * max_dist; }
  void set_max_sulfur_dist(double max_sulfur_dist) { max_sulfur_dist_sq = max_sulfur_dist * max_sulfur_dist; }
};

inline double get_atom_geometry_angle(HybridizationType hybridization) {
  switch (hybridization) {
  case HybridizationType::SP:
    return deg_to_rad(180.0);
  case HybridizationType::SP2:
    return deg_to_rad(120.0);
  case HybridizationType::SP3:
    return deg_to_rad(109.4721);
  case HybridizationType::SP3D2:
    return deg_to_rad(90.0);
  default:
    return deg_to_rad(120.0);
  }
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

AtomType add_hydrogen_donor(const RDKit::RWMol &mol, const RDKit::Atom &atom);
AtomType add_hydrogen_acceptor(const RDKit::RWMol &mol, const RDKit::Atom &atom);

// atom is in an aromatic ring with an electronegative element (N or O)
bool in_aromatic_ring_with_N_or_O(const RDKit::RWMol &mol, const RDKit::Atom &atom);

AtomType add_weak_hydrogen_donor(const RDKit::RWMol &mol, const RDKit::Atom &atom);

void find_hydrogen_bonds(Luni &luni, const GeometryOptions &opts, NSResults &neighbors, Contacts &container);

void find_weak_hydrogen_bonds(Luni &luni, const GeometryOptions &opts, NSResults &neighbors, Contacts &container);

} // namespace lahuta

#endif // LAHUTA_HYDROGEN_BONDS_HPP
