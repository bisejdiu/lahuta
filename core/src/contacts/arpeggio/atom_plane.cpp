#include <contacts/arpeggio/params.hpp>

#include "contacts/arpeggio/contacts.hpp"
#include "contacts/arpeggio/geo.hpp"
#include "elements.hpp"
#include "entities/context.hpp"

// clang-format off
namespace lahuta::arpeggio {

ContactRecipe<AtomRec, RingRec, DonorPiParams> make_donor_pi_recipe() {
  return {
    DonorPiParams{},
    +[](const AtomRec &rec) { return (rec.type & AtomType::HbondDonor) == AtomType::HbondDonor; },
    +[](const RingRec &rec) { return rec.aromatic; },
    +[](u32 a, u32 b, float d_sq, const ContactContext& ctx) -> InteractionType {

      const auto& params = ctx.get_params<DonorPiParams>();
      const auto rec_a = ctx.topology.atom(a);
      const auto rec_b = ctx.topology.ring(b);

      if (rec_a.atom.get().getIsAromatic()) return InteractionType::None; // FIX: ???

      auto angle = compute_angle(rec_b, ctx.conformer(), ctx.conformer().getAtomPos(rec_a.atom.get().getIdx()));
      if (!passes_angle_filter(angle, params.angle_cutoff)) return InteractionType::None;

      return InteractionType::DonorPi;
    }
  };
}

ContactRecipe<AtomRec, RingRec, SulphurPiParams> make_sulphur_pi_recipe() {
  return {
    SulphurPiParams{},
    +[](const AtomRec &rec) { return (rec.atom.get().getAtomicNum() == Element::S); },
    +[](const RingRec &rec) { return rec.aromatic; },
    +[](u32 a, u32 b, float d_sq, const ContactContext& ctx) -> InteractionType {

      const auto rec_a = ctx.topology.atom(a);
      const auto rec_b = ctx.topology.ring(b);

      auto res_a_info = static_cast<const RDKit::AtomPDBResidueInfo*>(rec_a.atom.get().getMonomerInfo());

      // NOTE: this can be added as part of the predicate
      if (!res_a_info || res_a_info->getResidueName() != "MET") return InteractionType::None;

      return InteractionType::SulphurPi;
    }
  };
}

ContactRecipe<AtomRec, RingRec, CarbonPiParams> make_carbon_pi_recipe() {
  return {
    CarbonPiParams{},
    +[](const AtomRec &rec) {
      bool is_C_atom = rec.atom.get().getAtomicNum() == Element::C;
      bool is_whd = (rec.type & AtomType::WeakHbondDonor) == AtomType::WeakHbondDonor;
      return is_C_atom && is_whd;
    },
    +[](const RingRec &rec) { return rec.aromatic; },
    +[](u32 a, u32 b, float d_sq, const ContactContext& ctx) -> InteractionType {

      const auto& params = ctx.get_params<CarbonPiParams>();
      const auto rec_a = ctx.topology.atom(a);
      const auto rec_b = ctx.topology.ring(b);

      if (rec_a.atom.get().getIsAromatic()) return InteractionType::None; // FIX: ???

      auto angle = compute_angle(rec_b, ctx.conformer(), ctx.conformer().getAtomPos(rec_a.atom.get().getIdx()));
      if (!passes_angle_filter(angle, params.angle_cutoff)) return InteractionType::None;

      // FIX: need to make sure we remove metals here
      return InteractionType::CarbonPi;
    }
  };
}

ContactRecipe<AtomRec, RingRec, CationPiParams> make_cation_pi_recipe() {
  return {
    CationPiParams{},
    +[](const AtomRec &rec) { return (rec.type & AtomType::PositiveCharge) == AtomType::PositiveCharge; },
    +[](const RingRec &rec) { return rec.aromatic; },
    +[](u32 a, u32 b, float d_sq, const ContactContext& ctx) -> InteractionType {

      const auto& params = ctx.get_params<CationPiParams>();
      const auto rec_a = ctx.topology.atom(a);
      const auto rec_b = ctx.topology.ring(b);

      if (rec_a.atom.get().getIsAromatic()) return InteractionType::None; // FIX: ???

      auto angle = compute_angle(rec_b, ctx.conformer(), ctx.conformer().getAtomPos(rec_a.atom.get().getIdx()));
      if (!passes_angle_filter(angle, params.angle_cutoff)) return InteractionType::None;

      return InteractionType::CationPi;
    }
  };
}

} // namespace lahuta::arpeggio
