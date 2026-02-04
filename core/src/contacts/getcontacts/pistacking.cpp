/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto concat = [](auto&&... args) {
 *     std::string result;
 *     ((result += std::string_view(args)), ...);
 *     return result;
 *   };
 *   return concat("besian", "sejdiu", "@gmail.com");
 * }();
 *
 */

#include <cmath>
#include <limits>
#include <optional>

#include "chemistry/utils.hpp"
#include "contacts/getcontacts/contacts.hpp"
#include "contacts/getcontacts/utils.hpp"
#include "entities/context.hpp"

// clang-format off
namespace lahuta::getcontacts {

namespace {

struct PiStackingEvaluation {
  double center_distance;
  double normal_angle;
  double psi_a;
  double psi_b;
};

std::optional<PiStackingEvaluation> evaluate_ring_pair(const ContactContext& ctx, const RingRec& ring_a, const RingRec& ring_b) {
  if (ring_a.atoms.empty() || ring_b.atoms.empty()) return std::nullopt;

  const auto& conf = ctx.conformer();

  const RDGeom::Point3D center_a = ring_a.center(conf);
  const RDGeom::Point3D center_b = ring_b.center(conf);

  const auto normal_a = ring_a.normal(conf);
  const auto normal_b = ring_b.normal(conf);

  const double dist = (center_a - center_b).length();
  if (dist <= std::numeric_limits<double>::min()) return std::nullopt;

  const double angle = detail::vector_angle_deg(normal_a, normal_b);
  const double psi_a = detail::psi_angle_deg(center_a, center_b, normal_a);
  const double psi_b = detail::psi_angle_deg(center_b, center_a, normal_b);

  return PiStackingEvaluation{dist, angle, psi_a, psi_b};
}

} // namespace

ContactRecipe<RingRec, RingRec, PiStackingParams> make_pi_stacking_recipe() {
  return {
    PiStackingParams{},
    +[](const RingRec& rec) { return rec.aromatic; },
    +[](const RingRec& rec) { return rec.aromatic; },
    +[](std::uint32_t idx_a, std::uint32_t idx_b, float /*dist_sq*/, const ContactContext& ctx) -> InteractionType {
      const auto& params = ctx.get_params<PiStackingParams>();
      const auto& ring_a = ctx.topology.ring(idx_a);
      const auto& ring_b = ctx.topology.ring(idx_b);

      if (is_same_residue(ctx.molecule(), ring_a.atoms.front(), ring_b.atoms.front())) return InteractionType::None;

      auto eval = evaluate_ring_pair(ctx, ring_a, ring_b);
      if (!eval) return InteractionType::None;

      if (eval->center_distance > params.centroid_cutoff) return InteractionType::None;

      // Parallel stacking: normals aligned (0 or 180) within cutoff
      const double aligned = std::min(std::fabs(eval->normal_angle), std::fabs(180.0 - eval->normal_angle));
      if (aligned > params.angle_cutoff_deg) return InteractionType::None;

      if (std::min(eval->psi_a, eval->psi_b) > params.psi_cutoff_deg) return InteractionType::None;
      return InteractionType::PiStackingP;
    }
  };
}

ContactRecipe<RingRec, RingRec, TStackingParams> make_t_stacking_recipe() {
  return {
    TStackingParams{},
    +[](const RingRec& rec) { return rec.aromatic; },
    +[](const RingRec& rec) { return rec.aromatic; },
    +[](std::uint32_t idx_a, std::uint32_t idx_b, float /*dist_sq*/, const ContactContext& ctx) -> InteractionType {
      const auto& params = ctx.get_params<TStackingParams>();
      const auto& ring_a = ctx.topology.ring(idx_a);
      const auto& ring_b = ctx.topology.ring(idx_b);

      if (is_same_residue(ctx.molecule(), ring_a.atoms.front(), ring_b.atoms.front())) return InteractionType::None;

      auto eval = evaluate_ring_pair(ctx, ring_a, ring_b);
      if (!eval) return InteractionType::None;

      if (eval->center_distance > params.centroid_cutoff) return InteractionType::None;

      const double deviation = std::fabs(eval->normal_angle - 90.0);
      if (deviation > params.angle_cutoff_deg) return InteractionType::None;

      if (std::min(eval->psi_a, eval->psi_b) > params.psi_cutoff_deg) return InteractionType::None;

      return InteractionType::PiStackingT;
    }
  };
}

} // namespace lahuta::getcontacts
