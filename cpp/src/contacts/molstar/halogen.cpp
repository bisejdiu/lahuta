#include "contacts/halo_geo_validity.hpp"
#include "contacts/molstar/contacts.hpp"

// clang-format off
namespace lahuta {

ContactRecipe<AtomRec, AtomRec, HalogenParams> make_halogen_recipe() {
  return {
    HalogenParams{},
    +[](const AtomRec &rec) { return (rec.type & AtomType::XbondDonor)    == AtomType::XbondDonor; },
    +[](const AtomRec &rec) { return (rec.type & AtomType::XBondAcceptor) == AtomType::XBondAcceptor; },
    // {params.distance_max, 0, 0, 0.7},
    +[](std::uint32_t rec_idx_a, std::uint32_t rec_idx_b, float dist, const ContactContext& ctx) -> InteractionType {
      const auto& params = ctx.get_params<HalogenParams>();
      const auto &donor    = ctx.topology.atom(rec_idx_a).atom;
      const auto &acceptor = ctx.topology.atom(rec_idx_b).atom;

      if (!halo_geo::are_geometrically_viable(ctx.molecule(), donor, acceptor, params)) return InteractionType::None;

      return InteractionType::Halogen;
    }
  };
}

} // namespace lahuta
