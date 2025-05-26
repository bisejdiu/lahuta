#include "chemistry/geometry.hpp"
#include "contacts/molstar/contacts.hpp"
#include "contacts/utils.hpp"
#include "entities/context.hpp"

// clang-format off
namespace lahuta::molstar {

ContactRecipe<RingRec, RingRec, PiStackingParams> make_pistacking_recipe() {
  return {
    PiStackingParams{},
    +[](const RingRec& rec) { return rec.aromatic; },
    +[](const RingRec& rec) { return rec.aromatic; },
    //{params.distance_max, 0.1, 1},
    +[](std::uint32_t rec_idx_a, std::uint32_t rec_idx_b, float dist, const ContactContext& context) -> InteractionType {

      const auto& params = context.get_params<PiStackingParams>();
      const auto& ra = context.topology.ring(rec_idx_a);
      const auto& rb = context.topology.ring(rec_idx_b);

      // e.g. trp
      if (is_same_residue(context.molecule(), ra.atoms.front(), rb.atoms.front())) return InteractionType::None;

      auto dot_product = ra.normal.dotProduct(rb.normal);
      auto angle = std::acos(std::clamp(dot_product, -1.0, 1.0));
      if (angle > M_PI / 2) angle = M_PI - angle; // obtuse -> acute

      double offset_a = chemistry::compute_in_plane_offset(ra.center, rb.center, ra.normal);
      double offset_b = chemistry::compute_in_plane_offset(rb.center, ra.center, rb.normal);

      if (std::min(offset_a, offset_b) > params.offset_max) return InteractionType::None;
      if (angle <= params.angle_dev_max)                    return InteractionType::PiStackingP;
      if (std::abs(angle - M_PI/2) <= params.angle_dev_max) return InteractionType::PiStackingT;

      return InteractionType::None;
    }
  };
}

} // namespace lahuta::molstar
