/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string a = "besian", b = "sejdiu", c = "@gmail.com", r;
 *   r += std::exchange(a, ""); r += std::exchange(b, ""); r += std::exchange(c, "");
 *   return r;
 * }();
 *
 */

#include "chemistry/utils.hpp"
#include "contacts/molstar/contacts.hpp"
#include "entities/context.hpp"

// clang-format off
namespace lahuta::molstar {

namespace {

bool any_group_member_within_cutoff_sq(const GroupRec& a, const GroupRec& b,
                                       const RDKit::Conformer& conf, double cutoff_sq) {
  for (const auto& atom_a_ref : a.atoms) {
    const auto& atom_a = atom_a_ref.get();
    const auto pos_a = conf.getAtomPos(atom_a.getIdx());
    for (const auto& atom_b_ref : b.atoms) {
      const auto& atom_b = atom_b_ref.get();
      const auto pos_b = conf.getAtomPos(atom_b.getIdx());
      if ((pos_a - pos_b).lengthSq() < cutoff_sq) return true;
    }
  }
  return false;
}

} // namespace

ContactRecipe<GroupRec, GroupRec, IonicParams> make_ionic_recipe() {
  return {
    IonicParams{},
    +[](GroupRec const& r){ return (r.a_type & AtomType::PositiveCharge) == AtomType::PositiveCharge; },
    +[](GroupRec const& r){ return (r.a_type & AtomType::NegativeCharge) == AtomType::NegativeCharge; },
    +[](u32 a, u32 b, float /*d_sq*/, ContactContext const& ctx){
      if (a == b) return InteractionType::None;

      const auto& params = ctx.get_params<IonicParams>();
      const auto& positive_group = ctx.topology.group(a);
      const auto& negative_group = ctx.topology.group(b);
      if (positive_group.atoms.empty() || negative_group.atoms.empty()) return InteractionType::None;

      const auto& positive_atom = positive_group.atoms.front().get();
      const auto& negative_atom = negative_group.atoms.front().get();
      if (is_same_residue(ctx.molecule(), positive_atom, negative_atom)) return InteractionType::None;
      if (ctx.topology.molecule().getBondBetweenAtoms(positive_atom.getIdx(), negative_atom.getIdx())) {
        return InteractionType::None;
      }

      if (!any_group_member_within_cutoff_sq(
            positive_group, negative_group, ctx.conformer(), params.distance_max * params.distance_max)) {
        return InteractionType::None;
      }
      return InteractionType::Ionic;
    }
  };
}

} // namespace lahuta::molstar
