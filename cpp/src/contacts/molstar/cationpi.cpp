#include "contacts/molstar/contacts.hpp"
#include "chemistry/geometry.hpp"
#include "contacts/utils.hpp"
#include "entities/context.hpp"

// clang-format off
namespace lahuta::molstar {

ContactRecipe<GroupRec, RingRec, CationPiParams> make_cationpi_recipe() {
  return {
    CationPiParams{},
    +[](const GroupRec& rec) { return (rec.a_type & AtomType::PositiveCharge) == AtomType::PositiveCharge; },
    +[](const RingRec & rec)  { return true; }, // we miss positive hits bc we miss genuine aromatic rings in our perception routine
    +[](std::uint32_t rec_idx_a, std::uint32_t rec_idx_b, float dist, const ContactContext& ctx) -> InteractionType {

      const auto& params = ctx.get_params<CationPiParams>();
      const auto &ci = ctx.topology.group(rec_idx_a);
      const auto &ri = ctx.topology.ring(rec_idx_b);

      if (is_same_residue(ctx.molecule(), ri.atoms.front(), ci.atoms.front())) return InteractionType::None;

      auto offset = chemistry::compute_in_plane_offset(ci.center, ri.center, ri.normal);
      if (offset > params.offset_max) return InteractionType::None;
      return InteractionType::CationPi;
    }
  };
}

} // namespace lahuta::molstar
