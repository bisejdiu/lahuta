#ifndef LAHUTA_CONTACTS_GEO_VALIDITY_HPP
#define LAHUTA_CONTACTS_GEO_VALIDITY_HPP

#include "chemistry/geometry.hpp"
#include "contacts/molstar/params.hpp"

#include "GraphMol/Atom.h"
#include "GraphMol/RWMol.h"

// clang-format off
namespace lahuta::molstar {


inline bool
are_geometrically_viable( const RDKit::RWMol &mol, const RDKit::Atom &donor, const RDKit::Atom &acceptor, const HalogenParams &opts) {

  auto [halogen_angles, _] = chemistry::calculate_angle(mol, donor, acceptor, true);
  if (halogen_angles.size() != 1) return false;

  if (opts.optimal_angle - halogen_angles[0] > opts.angle_max) return false;

  auto [acceptor_angles, __] = chemistry::calculate_angle(mol, acceptor, donor, true);
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

} // namespace lahuta::molstar

#endif // LAHUTA_CONTACTS_GEO_VALIDITY_HPP
