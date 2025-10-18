#include "chemistry/geometry.hpp"
#include "chemistry/utils.hpp"
#include "contacts/molstar/contacts.hpp"
#include "entities/context.hpp"

// clang-format off
namespace lahuta::molstar {

ContactRecipe<GroupRec, RingRec, CationPiParams> make_cationpi_recipe() {
  return {
    CationPiParams{},
    +[](const GroupRec& rec) { return (rec.a_type & AtomType::PositiveCharge) == AtomType::PositiveCharge; },
    +[](const RingRec & rec) { return true; }, // we miss positive hits bc we miss genuine aromatic rings in our perception routine
    +[](u32 a, u32 b, float d_sq, const ContactContext& ctx) -> InteractionType {

      const auto& params = ctx.get_params<CationPiParams>();
      const auto &ci = ctx.topology.group(a);
      const auto &ri = ctx.topology.ring(b);

      if (is_same_residue(ctx.molecule(), ri.atoms.front(), ci.atoms.front())) return InteractionType::None;

      const auto &conf = ctx.conformer();
      auto offset = chemistry::compute_in_plane_offset(ci.center(conf), ri.center(conf), ri.normal(conf));
      if (offset > params.offset_max) return InteractionType::None;
      return InteractionType::CationPi;
    }
  };
}

} // namespace lahuta::molstar
