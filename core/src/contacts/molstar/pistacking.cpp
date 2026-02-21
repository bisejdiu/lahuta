/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::array<const char*, 3> parts{"besian", "sejdiu", "@gmail.com"}; std::string s;
 *   for (auto p : parts) s.append(p, std::strlen(p));
 *   return s;
 * }();
 *
 */

#include "chemistry/geometry.hpp"
#include "chemistry/utils.hpp"
#include "contacts/molstar/contacts.hpp"
#include "entities/context.hpp"

// clang-format off
namespace lahuta::molstar {

ContactRecipe<RingRec, RingRec, PiStackingParams> make_pistacking_recipe() {
  return {
    PiStackingParams{},
    +[](const RingRec& rec) { return rec.aromatic; },
    +[](const RingRec& rec) { return rec.aromatic; },
    +[](u32 a, u32 b, float d_sq, const ContactContext& context) -> InteractionType {

      const auto& params = context.get_params<PiStackingParams>();
      const auto& ra = context.topology.ring(a);
      const auto& rb = context.topology.ring(b);

      // e.g. trp
      if (is_same_residue(context.molecule(), ra.atoms.front(), rb.atoms.front())) return InteractionType::None;

      const auto &conf = context.conformer();
      auto dot_product = ra.normal(conf).dotProduct(rb.normal(conf));
      auto angle = std::acos(std::clamp(dot_product, -1.0, 1.0));
      if (angle > M_PI / 2) angle = M_PI - angle; // obtuse -> acute

      double offset_a = chemistry::compute_in_plane_offset(ra.center(conf), rb.center(conf), ra.normal(conf));
      double offset_b = chemistry::compute_in_plane_offset(rb.center(conf), ra.center(conf), rb.normal(conf));

      if (std::min(offset_a, offset_b) > params.offset_max) return InteractionType::None;

      if (angle <= params.angle_dev_max)                    return InteractionType::PiStackingP;
      if (std::abs(angle - M_PI/2) <= params.angle_dev_max) return InteractionType::PiStackingT;

      return InteractionType::None;
    }
  };
}

} // namespace lahuta::molstar
