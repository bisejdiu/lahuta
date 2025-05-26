#include "contacts/arpeggio/contacts.hpp"
#include "contacts/utils.hpp"
#include "entities/context.hpp"
#include <contacts/arpeggio/params.hpp>

// clang-format off
namespace lahuta::arpeggio {

ContactRecipe<AtomRec, AtomRec, AromaticParams> make_aromatic_recipe() {
  return {
    AromaticParams{},
    +[](const AtomRec& rec) { return (rec.type & AtomType::Aromatic) == AtomType::Aromatic; },
    +[](std::uint32_t rec_idx_a, std::uint32_t rec_idx_b, float dist, const ContactContext& ctx) -> InteractionType {

      const auto& ring_a = ctx.topology.atom(rec_idx_a);
      const auto& ring_b = ctx.topology.atom(rec_idx_b);

      if (are_residueids_close(ctx.molecule(), ring_a.atom.get(), ring_b.atom.get(), 1)) return InteractionType::None;
      return InteractionType::PiStackingP;
    }
  };
}

} // namespace lahuta::arpeggio
