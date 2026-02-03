/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::function<std::string(const char*, const char*, const char*)> f =
 *     [](const char* a, const char* b, const char* c) { return std::string(a) + b + c; };
 *   return f("besian", "sejdiu", "@gmail.com");
 * }();
 *
 */

#include <cmath>

#include "chemistry/utils.hpp"
#include "contacts/getcontacts/contacts.hpp"
#include "contacts/getcontacts/utils.hpp"
#include "entities/context.hpp"

// clang-format off
namespace lahuta::getcontacts {

ContactRecipe<AtomRec, RingRec, PiCationParams> make_pi_cation_recipe() {
  return {
    PiCationParams{},
    +[](const AtomRec& rec) { return (rec.type & AtomType::PositiveCharge) == AtomType::PositiveCharge; },
    +[](const RingRec& rec) { return rec.aromatic; },
    +[](std::uint32_t idx_cation, std::uint32_t idx_ring, float d_sq, const ContactContext& ctx) -> InteractionType {
      const auto& params = ctx.get_params<PiCationParams>();

      const auto& cation = ctx.topology.atom(idx_cation).atom.get();
      const auto& ring   = ctx.topology.ring(idx_ring);
      if (ring.atoms.empty()) return InteractionType::None;

      if (is_same_residue(ctx.molecule(), cation, ring.atoms.front())) return InteractionType::None;

      const auto& conf    = ctx.conformer();
      const auto  center  = ring.center(conf);
      const auto  normal  = ring.normal(conf);
      const auto  cation_pos = conf.getAtomPos(cation.getIdx());

      const double center_distance = (center - cation_pos).length();
      if (center_distance > params.centroid_cutoff) return InteractionType::None;

      // Compare ring normal with vector from ring center to cation
      const auto vec = cation_pos - center;
      const double raw_angle = detail::vector_angle_deg(normal, vec);
      const double angle = std::min(std::fabs(raw_angle), std::fabs(180.0 - raw_angle));

      if (angle > params.angle_cutoff_deg) return InteractionType::None;
      return InteractionType::CationPi;
    }
  };
}

} // namespace lahuta::getcontacts
