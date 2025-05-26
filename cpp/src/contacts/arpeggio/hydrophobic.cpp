#include "contacts/arpeggio/contacts.hpp"
#include "contacts/utils.hpp"
#include "entities/context.hpp"
#include <contacts/arpeggio/params.hpp>

// clang-format off
namespace lahuta::arpeggio {

ContactRecipe<AtomRec, AtomRec, HydrophobicParams> make_hydrophobic_recipe() {
  return {
    HydrophobicParams{},
    +[](const AtomRec &rec) { return (rec.type & AtomType::Hydrophobic) == AtomType::Hydrophobic; },
    +[](std::uint32_t rec_idx_a, std::uint32_t rec_idx_b, float dist, const ContactContext& ctx) -> InteractionType {

      const auto& params = ctx.get_params<HydrophobicParams>();
      const auto rec_a = ctx.topology.atom(rec_idx_a);
      const auto rec_b = ctx.topology.atom(rec_idx_b);

      if (are_residueids_close(ctx.molecule(), rec_a.atom, rec_b.atom, 1)) return InteractionType::None;
      return InteractionType::Hydrophobic;
    }
  };
}

ContactRecipe<AtomRec, AtomRec, VanDerWaalsParams> make_vdw_recipe() {
  return {
    VanDerWaalsParams{},
    +[](const AtomRec& rec) { return true; },
    +[](const AtomRec& rec) { return true; },
    +[](std::uint32_t rec_idx_a, std::uint32_t rec_idx_b, float dist_sq, const ContactContext& ctx) -> InteractionType {

      const auto& params = ctx.get_params<VanDerWaalsParams>();
      const auto& atom_a = ctx.topology.atom(rec_idx_a).atom.get();
      const auto& atom_b = ctx.topology.atom(rec_idx_b).atom.get();

      if (atom_a.getAtomicNum() == Element::H || atom_b.getAtomicNum() == Element::H) return InteractionType::None;

      auto vdw_a = vdw_radius(static_cast<Element>(atom_a.getAtomicNum()));
      auto vdw_b = vdw_radius(static_cast<Element>(atom_b.getAtomicNum()));

      // FIX: remove clashes is not used.
      float sum_vdw = vdw_a + vdw_b;
      float max_sq = (sum_vdw + params.vdw_comp_factor) * (sum_vdw + params.vdw_comp_factor);

      if (dist_sq > max_sq || dist_sq < (sum_vdw * sum_vdw))  return InteractionType::None;
      if (are_residueids_close(ctx.molecule(), atom_a, atom_b, 1)) return InteractionType::None;

      return InteractionType::VanDerWaals;
    }
  };
}

} // namespace lahuta::arpeggio
