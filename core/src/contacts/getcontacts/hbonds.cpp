/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::array<std::string_view, 3> parts{"besian", "sejdiu", "@gmail.com"}; std::string s;
 *   for (auto p : parts) s.append(p);
 *   return s;
 * }();
 *
 */

#include <cmath>
#include <optional>

#include <rdkit/GraphMol/RWMol.h>

#include "contacts/getcontacts/contacts.hpp"
#include "contacts/getcontacts/utils.hpp"
#include "chemistry/elements.hpp"
#include "entities/context.hpp"

// clang-format off
namespace lahuta::getcontacts {

namespace {

std::optional<double> best_dha_angle_deg(const ContactContext& ctx, const RDKit::Atom& donor, const RDKit::Atom& acceptor) {
  const auto& mol  = ctx.molecule();
  const auto& conf = ctx.conformer();

  double best_angle = -1.0;
  bool   found      = false;

  for (const auto* nbr : mol.atomNeighbors(&donor)) {
    if (!nbr) continue;
    if (nbr->getIdx() == acceptor.getIdx()) continue;
    if (nbr->getAtomicNum() != Element::H) continue;

    const auto& h_pos  = conf.getAtomPos(nbr->getIdx());
    const auto& d_pos  = conf.getAtomPos(donor.getIdx());
    const auto& a_pos  = conf.getAtomPos(acceptor.getIdx());

    RDGeom::Point3D vec_d = d_pos - h_pos;
    RDGeom::Point3D vec_a = a_pos - h_pos;
    const double angle = detail::vector_angle_deg(vec_d, vec_a);
    if (angle > best_angle) {
      best_angle = angle;
      found = true;
    }
  }

  if (!found) return std::nullopt;
  return best_angle;
}

} // namespace

ContactRecipe<AtomRec, AtomRec, HydrogenBondParams> make_hbond_recipe() {
  return {
    HydrogenBondParams{},
    +[](const AtomRec& rec) { return (rec.type & AtomType::HbondDonor)    == AtomType::HbondDonor; },
    +[](const AtomRec& rec) { return (rec.type & AtomType::HbondAcceptor) == AtomType::HbondAcceptor; },
    +[](std::uint32_t idx_d, std::uint32_t idx_a, float d_sq, const ContactContext& ctx) -> InteractionType {

      const auto& params   = ctx.get_params<HydrogenBondParams>();
      const auto& donor    = ctx.topology.atom(idx_d).atom.get();
      const auto& acceptor = ctx.topology.atom(idx_a).atom.get();

      // Exclude S-involving H-bonds if parameter is set
      if (params.exclude_sulfur_hbonds) {
        if (donor.getAtomicNum() == Element::S || acceptor.getAtomicNum() == Element::S) {
          return InteractionType::None;
        }
      }

      // Allow water-mediated interactions
      // if (!params.include_water && detail::is_water(donor) && detail::is_water(acceptor)) return InteractionType::None;

      if (detail::residues_too_close(ctx, donor, acceptor, params.min_residue_offset)) return InteractionType::None;

      const double dist = std::sqrt(d_sq);
      if (dist > params.distance_cutoff) return InteractionType::None;

      auto angle_opt = best_dha_angle_deg(ctx, donor, acceptor);
      if (angle_opt) {
        const double angle = *angle_opt;
        const double min_angle = 180.0 - params.angle_tolerance_deg;
        if (angle < min_angle) return InteractionType::None;
      }

      return InteractionType::HydrogenBond;
    }
  };
}

} // namespace lahuta::getcontacts
