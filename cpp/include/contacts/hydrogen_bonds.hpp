/*
    This code adapts functionality from the Mol* (mol-star) molecular visualization package,
    originally in TypeScript, for use in C++ with RDKit in the Lahuta project.

    Mol* is licensed under the MIT License (https://github.com/molstar/molstar),
    while this adaptation is released under the GNU General Public License (GPL).
*/
#ifndef LAHUTA_HYDROGEN_BONDS_HPP
#define LAHUTA_HYDROGEN_BONDS_HPP

#include "typing/types.hpp"
#include "entities/contact.hpp"
#include <GraphMol/Atom.h>
#include <GraphMol/MonomerInfo.h>
#include <GraphMol/RWMol.h>

// clang-format off
namespace lahuta {
class Topology;

constexpr double MAX_DIST = 3.5;
constexpr double MAX_SULFUR_DIST = 4.1;
constexpr double MAX_ACC_ANGLE_DEV = M_PI / 4.0;          // 45 degrees
constexpr double MAX_DON_ANGLE_DEV = M_PI / 4.0;          // 45 degrees
constexpr double MAX_DON_OUT_OF_PLANE_ANGLE = M_PI / 4.0; // 45 degrees
constexpr double MAX_ACC_OUT_OF_PLANE_ANGLE = M_PI / 2.0; // 90 degrees

struct HBondParameters {
  bool   ignore_hydrogens  = false;             // Ignore hydrogens in geometry calculations
  bool   include_backbone  = true;              // Include backbone-to-backbone hydrogen bonds
  bool   include_water     = true;              // Include water-to-water hydrogen bonds
  double max_dist          = MAX_DIST;          // Maximum distance for hydrogen bonds
  double max_sulfur_dist   = MAX_SULFUR_DIST;   // Maximum distance for sulfur atoms
  double max_don_angle_dev = MAX_DON_ANGLE_DEV; // Maximum deviation from ideal donor angle
  double max_acc_angle_dev = MAX_ACC_ANGLE_DEV; // Maximum deviation from ideal acceptor angle
  double max_don_out_of_plane_angle = MAX_DON_OUT_OF_PLANE_ANGLE; // Maximum out-of-plane deviation for donor
  double max_acc_out_of_plane_angle = MAX_ACC_OUT_OF_PLANE_ANGLE; // Maximum out-of-plane deviation for acceptor

};

AtomType add_hydrogen_donor     (const RDKit::RWMol &mol, const RDKit::Atom &atom);
AtomType add_hydrogen_acceptor  (const RDKit::RWMol &mol, const RDKit::Atom &atom);
AtomType add_weak_hydrogen_donor(const RDKit::RWMol &mol, const RDKit::Atom &atom);

class Luni;
ContactSet find_hydrogen_bonds     (const Topology &topology, const HBondParameters &opts = HBondParameters{});
ContactSet find_weak_hydrogen_bonds(const Topology &topology, HBondParameters opts = HBondParameters{});

} // namespace lahuta

#endif // LAHUTA_HYDROGEN_BONDS_HPP
