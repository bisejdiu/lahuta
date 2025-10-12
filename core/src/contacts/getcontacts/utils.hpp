#ifndef LAHUTA_CONTACTS_GETCONTACTS_UTILS_HPP
#define LAHUTA_CONTACTS_GETCONTACTS_UTILS_HPP

#include <algorithm>
#include <cmath>
#include <limits>

#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/MonomerInfo.h>
#include <GraphMol/RWMol.h>

#include "chemistry/utils.hpp"
#include "definitions.hpp"
#include "elements.hpp"
#include "entities/context.hpp"

// clang-format off
namespace lahuta::getcontacts::detail {


inline double rad_to_deg(double rad) noexcept {
  constexpr double Pi = 3.14159265358979323846264338327950288;
  return rad * (180.0 / Pi);
}

inline double clamp_cosine(double value) noexcept {
  if (value >  1.0) return  1.0;
  if (value < -1.0) return -1.0;
  return value;
}

inline double vector_angle_deg(const RDGeom::Point3D& v1, const RDGeom::Point3D& v2) noexcept {
  const auto len1 = v1.length();
  const auto len2 = v2.length();
  if (len1 <= std::numeric_limits<double>::min() || len2 <= std::numeric_limits<double>::min()) return 180.0;
  const double cos_val = clamp_cosine(v1.dotProduct(v2) / (len1 * len2));
  return rad_to_deg(std::acos(cos_val));
}

inline double psi_angle_deg(const RDGeom::Point3D& center_a, const RDGeom::Point3D& center_b, const RDGeom::Point3D& normal_a) noexcept {
  const RDGeom::Point3D vec = center_b - center_a;
  const double angle = vector_angle_deg(normal_a, vec);
  return std::min(std::fabs(angle), std::fabs(180.0 - angle));
}

// Check if two atoms are directly bonded sulfur atoms (S-S disulfide bond)
inline bool is_disulfide_pair(const RDKit::RWMol& mol, const RDKit::Atom& a, const RDKit::Atom& b) noexcept {
  if (a.getAtomicNum() != Element::S || b.getAtomicNum() != Element::S) return false;
  return mol.getBondBetweenAtoms(a.getIdx(), b.getIdx()) != nullptr;
}

// Check if two atoms belong to CYS residues that form a disulfide bridge.
inline bool is_cys_disulfide_contact(const RDKit::RWMol& mol, const RDKit::Atom& a, const RDKit::Atom& b) noexcept {
  // Check if both atoms are from CYS residues
  auto info_a = static_cast<const RDKit::AtomPDBResidueInfo*>(a.getMonomerInfo());
  auto info_b = static_cast<const RDKit::AtomPDBResidueInfo*>(b.getMonomerInfo());

  if (!info_a || !info_b) return false;
  if (info_a->getResidueName() != "CYS" || info_b->getResidueName() != "CYS") return false;

  // If they're from the same residue, it's not a disulfide contact
  if (info_a->getResidueNumber() == info_b->getResidueNumber() && info_a->getChainId() == info_b->getChainId()) return false;

  const RDKit::Atom* sg_a = nullptr;
  const RDKit::Atom* sg_b = nullptr;

  for (const auto& atom : mol.atoms()) {
    auto info = static_cast<const RDKit::AtomPDBResidueInfo*>(atom->getMonomerInfo());
    if (!info) continue;

    if (info->getResidueNumber() == info_a->getResidueNumber() &&
        info->getChainId() == info_a->getChainId() &&
        info->getName() == "SG") {
      sg_a = atom;
    }

    if (info->getResidueNumber() == info_b->getResidueNumber() &&
        info->getChainId() == info_b->getChainId() &&
        info->getName() == "SG") {
      sg_b = atom;
    }
  }

  if (sg_a && sg_b) {
    return mol.getBondBetweenAtoms(sg_a->getIdx(), sg_b->getIdx()) != nullptr;
  }

  return false;
}

inline bool is_water(const RDKit::Atom& atom) noexcept {
  auto info = static_cast<const RDKit::AtomPDBResidueInfo*>(atom.getMonomerInfo());
  if (!info) return false;
  const auto resn = info->getResidueName();
  return resn == "HOH" || resn == "H2O" || resn == "WAT" || resn == "SOL" || resn == "TIP3" || resn == "TIP";
}

inline bool is_ligand(const RDKit::Atom& atom) noexcept {
  auto info = static_cast<const RDKit::AtomPDBResidueInfo*>(atom.getMonomerInfo());
  if (!info) return true;  // If no residue info, assume it's a ligand
  const auto resn = info->getResidueName();
  return !definitions::is_standard_protein(resn) &&
         !definitions::is_base(resn) &&
         !definitions::is_polymer(resn) &&
         !is_water(atom);
}

inline bool is_hydrophobic_carbon(const RDKit::RWMol& mol, const RDKit::Atom& atom) noexcept {
  if (atom.getAtomicNum() != Element::C) return false;
  for (const auto& bond : mol.atomBonds(&atom)) {
    const RDKit::Atom* nbr = bond->getOtherAtom(&atom);
    const auto z = nbr->getAtomicNum();
    if (z != Element::C && z != Element::H) return false;
  }
  return true;
}

inline bool residues_too_close(const ContactContext& ctx, const RDKit::Atom& a, const RDKit::Atom& b, int min_offset) noexcept {
  if (min_offset <= 0) return false;
  return are_residueids_close(ctx.molecule(), a, b, min_offset - 1);
}

} // namespace lahuta::getcontacts::detail

#endif // LAHUTA_CONTACTS_GETCONTACTS_UTILS_HPP
