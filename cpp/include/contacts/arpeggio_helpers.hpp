#ifndef LAHUTA_ARPEGGIO_HELPERS_HPP
#define LAHUTA_ARPEGGIO_HELPERS_HPP

#include "GraphMol/Atom.h"
#include "elements.hpp"
#include <entities/records.hpp>

// clang-format off
namespace lahuta {

// Find hydrogen atoms bonded to a donor atom
inline std::vector<const RDKit::Atom*> find_bonded_hydrogens(const RDKit::ROMol& mol, const RDKit::Atom* donor_atom) {
  std::vector<const RDKit::Atom*> hydrogens;

  for (const auto& bond : mol.atomBonds(donor_atom)) {
    const RDKit::Atom* neighbor = bond->getOtherAtom(donor_atom);
    if (neighbor->getAtomicNum() == Element::H) {
        hydrogens.push_back(neighbor);
    }
  }

  return hydrogens;
}

inline bool passes_hbond_distance_filter(
  const RDKit::ROMol& mol,
  const RDKit::Conformer& conformer,
  const RDKit::Atom* donor_atom,
  const RDKit::Atom* acceptor_atom,
  double vdw_comp_factor = 0.1) {
    auto donor_hydrogens = find_bonded_hydrogens(mol, donor_atom);
    if (donor_hydrogens.empty()) return false;

    auto acceptor_vdw = vdw_radius(static_cast<Element>(acceptor_atom->getAtomicNum()));
    auto hydrogen_vdw = vdw_radius(Element::H);
    double max_distance = acceptor_vdw + hydrogen_vdw + vdw_comp_factor;

    RDGeom::Point3D acceptor_pos = conformer.getAtomPos(acceptor_atom->getIdx());

    for (const auto* hydrogen : donor_hydrogens) {
      RDGeom::Point3D hydrogen_pos = conformer.getAtomPos(hydrogen->getIdx());

      double distance = (hydrogen_pos - acceptor_pos).length(); // FIX: we compute sqrt
      if (distance <= max_distance) return true;
    }

    return false;
}

// Compute angle at "vertex" between points p1 and p2
inline double compute_angle_rad(const RDGeom::Point3D &vertex, const RDGeom::Point3D &p1, const RDGeom::Point3D &p2) {
  RDGeom::Point3D v1 = p1 - vertex;
  RDGeom::Point3D v2 = p2 - vertex;
  v1.normalize();
  v2.normalize();
  double dot = v1.dotProduct(v2);
  dot = std::clamp(dot, -1.0, 1.0);
  return std::acos(dot);
}


inline bool passes_hbond_angle_filter(
  const RDKit::ROMol& mol,
  const RDKit::Conformer& conformer,
  const RDKit::Atom* donor_atom,
  const RDKit::Atom* acceptor_atom,
  double angle_cutoff_rad,
  bool weak = false) {
    auto donor_hydrogens = find_bonded_hydrogens(mol, donor_atom);
    if (donor_hydrogens.empty()) return false;

    RDGeom::Point3D donor_pos = conformer.getAtomPos(donor_atom->getIdx());
    RDGeom::Point3D acceptor_pos = conformer.getAtomPos(acceptor_atom->getIdx());

    for (const auto* hydrogen : donor_hydrogens) {
      RDGeom::Point3D hydrogen_pos = conformer.getAtomPos(hydrogen->getIdx());

      // Calculate D-H...A angle (at H)
      double angle = compute_angle_rad(donor_pos, hydrogen_pos, acceptor_pos);
      if (angle >= angle_cutoff_rad) return true;
    }

    return false;
}


inline double compute_angle(const RingRec &rd, const RDGeom::Point3D &point) {
  auto vector_point_to_plane = point - rd.center;
  vector_point_to_plane.normalize();

  double cos_theta = vector_point_to_plane.dotProduct(rd.normal);
  cos_theta = std::max(-1.0, std::min(1.0, cos_theta));

  double theta_radians = std::acos(cos_theta); // in radians
  return theta_radians * (180.0 / M_PI);
}

inline bool passes_angle_filter(double angle, double cutoff) {
  if (angle > 90.0) angle = 180.0 - angle; // “fold” angles >90 back in [0,90]
  return angle <= cutoff;
}

} // namespace lahuta

#endif // LAHUTA_ARPEGGIO_HELPERS_HPP
