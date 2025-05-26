#include "contacts/molstar/contacts.hpp"
#include "entities/context.hpp"
#include "typing/flags.hpp"

// clang-format off
namespace lahuta::molstar {

bool is_metalic(AtomType at1, AtomType at2) {
  using AtomTypeFlags::has;
  if (has(at1, AtomType::TransitionMetal)) return has(at2, AtomType::DativeBondPartner);
  if (has(at1, AtomType::IonicTypeMetal))  return has(at2, AtomType::IonicTypePartner);
  return false;
}

ContactRecipe<AtomRec, AtomRec, MetalicParams> make_metalic_recipe() {
  return {
    MetalicParams{},
    +[](const AtomRec &rec) { return (rec.type & (AtomType::IonicTypeMetal   | AtomType::TransitionMetal))   != AtomType::None; },
    +[](const AtomRec &rec) { return (rec.type & (AtomType::IonicTypePartner | AtomType::DativeBondPartner)) != AtomType::None; },
    +[](u32 a, u32 b, float d_sq, const ContactContext& ctx) -> InteractionType {
      const auto& opts = ctx.get_params<MetalicParams>();
      const auto &m  = ctx.topology.atom(a);
      const auto &mb = ctx.topology.atom(b);

      if (d_sq < opts.distance_max) return InteractionType::None;
      if (a == b) return InteractionType::None;

      if (!is_metalic(m.type, mb.type) && !is_metalic(mb.type, m.type)) return InteractionType::None;
      if (ctx.topology.molecule().getBondBetweenAtoms(m.atom.get().getIdx(), mb.atom.get().getIdx())) return InteractionType::None;

      return InteractionType::MetalCoordination;
    }
  };
}

} // namespace lahuta::molstar
