#include <cmath>

#include "contacts/getcontacts/contacts.hpp"
#include "contacts/getcontacts/utils.hpp"
#include "elements.hpp"
#include "entities/context.hpp"

// clang-format off
namespace lahuta::getcontacts {

ContactRecipe<AtomRec, AtomRec, VanDerWaalsParams> make_vdw_recipe() {
  return {
    VanDerWaalsParams{},
    +[](const AtomRec& rec) { return rec.atom.get().getAtomicNum() != Element::H; },
    +[](const AtomRec& rec) { return rec.atom.get().getAtomicNum() != Element::H; },
    +[](std::uint32_t idx_a, std::uint32_t idx_b, float dist_sq, const ContactContext& ctx) -> InteractionType {
      const auto& params = ctx.get_params<VanDerWaalsParams>();
      const auto& mol    = ctx.molecule();

      const auto& atom_a = ctx.topology.atom(idx_a).atom.get();
      const auto& atom_b = ctx.topology.atom(idx_b).atom.get();

      if (detail::residues_too_close(ctx, atom_a, atom_b, params.min_residue_offset)) return InteractionType::None;
      if (detail::is_cys_disulfide_contact(mol, atom_a, atom_b)) return InteractionType::None;

      const double distance = std::sqrt(dist_sq);
      const double vdw_a    = elements::vdw_radius(static_cast<Element>(atom_a.getAtomicNum()));
      const double vdw_b    = elements::vdw_radius(static_cast<Element>(atom_b.getAtomicNum()));
      const double cutoff   = vdw_a + vdw_b + params.epsilon;
      if (distance > cutoff) return InteractionType::None;

      return InteractionType::VanDerWaals;
    }
  };
}

} // namespace lahuta::getcontacts
