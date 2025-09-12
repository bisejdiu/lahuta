#include "contacts/arpeggio/contacts.hpp"
#include "chemistry/utils.hpp"
#include "entities/context.hpp"
#include <contacts/arpeggio/params.hpp>

// clang-format off
namespace lahuta::arpeggio {

ContactRecipe<AtomRec, AtomRec, CarbonylParams> make_carbonyl_recipe() {
  return {
    CarbonylParams{},
    +[](const AtomRec& rec) { return (rec.type & AtomType::CarbonylCarbon) == AtomType::CarbonylCarbon; },
    +[](const AtomRec& rec) { return (rec.type & AtomType::CarbonylOxygen) == AtomType::CarbonylOxygen; },
    +[](u32 a, u32 b, float d_sq, const ContactContext& ctx) -> InteractionType {

      const auto& rec_a = ctx.topology.atom(a);
      const auto& rec_b = ctx.topology.atom(b);

      if (are_residueids_close(ctx.molecule(), rec_a.atom, rec_b.atom, 1)) return InteractionType::None;
      return InteractionType::Carbonyl;
    }
  };
}

} // namespace lahuta::arpeggio
