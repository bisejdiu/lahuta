#ifndef LAHUTA_CONTACTS_HBOND_GEO_VALIDITY_HPP
#define LAHUTA_CONTACTS_HBOND_GEO_VALIDITY_HPP

#include "chemistry/geometry.hpp"
#include "chemistry/neighbors.hpp"
#include "contacts/molstar/contacts.hpp"

// clang-format off
namespace lahuta::hb_geo {

inline bool
are_geometrically_viable(const RDKit::RWMol &mol, const RDKit::Atom &donor, const RDKit::Atom &acceptor, const HBondParameters &opts) {

  // donor angles
  const auto &[don_angles, don_h_angles] = chemistry::calculate_angle(mol, donor, acceptor, opts.ignore_hydrogens);
  double ideal_don_angle = chemistry::get_atom_geometry_angle(donor.getHybridization());

  auto exceeds_dev = [ideal_don_angle, &opts](double angle) { return std::abs(ideal_don_angle - angle) > opts.max_don_angle_dev; };
  if (std::any_of(don_angles.begin(), don_angles.end(), exceeds_dev)) return false;

  auto is_within_dev = [&opts](double h_angle) { return h_angle >= opts.max_don_angle_dev; };
  if (!don_h_angles.empty() && std::all_of(don_h_angles.begin(), don_h_angles.end(), is_within_dev)) return false;

  // out-of-plane angle for sp2 hybridized atoms
  if (donor.getHybridization() == HybridizationType::SP2) {
    auto out_of_plane = chemistry::compute_plane_angle(mol, donor, acceptor);
    if (out_of_plane.value_or(-1.0) > opts.max_don_out_of_plane_angle) return false;
  }

  // If hydrogens are being considered and any exist, use the nearest H, otherwise use donor
  const auto &donor_atom_for_acc = (!opts.ignore_hydrogens && !don_h_angles.empty())
          ? chemistry::find_closest_hydrogen_atom(mol, donor, acceptor).value_or(std::cref(donor)).get()
          : donor;

  // acceptor angles
  const auto &[acc_angles, acc_h_angles] = chemistry::calculate_angle(mol, acceptor, donor_atom_for_acc, opts.ignore_hydrogens);
  double ideal_acc_angle = chemistry::get_atom_geometry_angle(acceptor.getHybridization());

  // Check acceptor angles (do not limit large acceptor angles)
  double min_acc = ideal_acc_angle - opts.max_acc_angle_dev;
  auto is_too_large = [min_acc](double angle) { return angle < min_acc; };
  if (std::any_of(acc_angles  .begin(), acc_angles  .end(), is_too_large) ||
      std::any_of(acc_h_angles.begin(), acc_h_angles.end(), is_too_large)) {
    return false;
  }

  // out-of-plane angle for sp2 hybridized atoms
  if (acceptor.getHybridization() == RDKit::Atom::HybridizationType::SP2) {
    auto out_of_plane = chemistry::compute_plane_angle(mol, acceptor, donor_atom_for_acc);
    if (out_of_plane.value_or(-1.0) > opts.max_acc_out_of_plane_angle) return false;
  }
  return true;
}

} // namespace lahuta::contacts

#endif // LAHUTA_CONTACTS_HBOND_GEO_VALIDITY_HPP
