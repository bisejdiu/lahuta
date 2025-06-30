#include "contacts/molstar/contacts.hpp"
#include "contacts/molstar/halo_geo_validity.hpp"
#include "entities/context.hpp"

// clang-format off
namespace lahuta::molstar {

ContactRecipe<AtomRec, AtomRec, HalogenParams> make_halogen_recipe() {
  return {
    HalogenParams{},
    +[](const AtomRec &rec) { return (rec.type & AtomType::XbondDonor)    == AtomType::XbondDonor; },
    +[](const AtomRec &rec) { return (rec.type & AtomType::XBondAcceptor) == AtomType::XBondAcceptor; },
    +[](u32 a, u32 b, float d_sq, const ContactContext& ctx) -> InteractionType {
      const auto& params = ctx.get_params<HalogenParams>();
      const auto &donor    = ctx.topology.atom(a).atom;
      const auto &acceptor = ctx.topology.atom(b).atom;

      if (!are_geometrically_viable(ctx.molecule(), donor, acceptor, params)) return InteractionType::None;

      return InteractionType::Halogen;
    }
  };
}

} // namespace lahuta::molstar
