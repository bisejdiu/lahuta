/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   constexpr std::string_view first = "besian", last = "sejdiu", host = "gmail.com";
 *   return std::string(first) + std::string(last) + "@" + std::string(host);
 * }();
 *
 */

#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/MonomerInfo.h>
#include <rdkit/GraphMol/RWMol.h>
#include "chemistry/common.hpp"
#include "chemistry/utils.hpp"
#include "contacts/molstar/contacts.hpp"
#include "contacts/molstar/hbond_geo_validity.hpp"
#include "chemistry/elements.hpp"
#include "entities/context.hpp"
#include "residues/definitions.hpp"

// clang-format off
namespace lahuta::molstar {

const std::array<std::string, 11> WaterResidues = {"HOH", "W", "SOL", "TIP3", "SPC", "H2O", "TIP4", "TIP", "DOD", "D3O", "WAT"};

bool is_water(const RDKit::Atom &atom) {
  auto res_info = static_cast<const RDKit::AtomPDBResidueInfo *>(atom.getMonomerInfo());
  if (!res_info) return false;
  return common::contains(WaterResidues, res_info->getResidueName());
}

bool is_water_hbond(const RDKit::Atom &atom_a, const RDKit::Atom &atom_b) {
  return is_water(atom_a) && is_water(atom_b);
}

bool are_molstar_residueids_close(
    const RDKit::ROMol &, const RDKit::Atom &atom_a, const RDKit::Atom &atom_b, int threshold) {
  auto info_a = static_cast<const RDKit::AtomPDBResidueInfo *>(atom_a.getMonomerInfo());
  auto info_b = static_cast<const RDKit::AtomPDBResidueInfo *>(atom_b.getMonomerInfo());

  if (!info_a || !info_b) return false;
  if (info_a->getChainId() != info_b->getChainId()) return false;
  if (info_a->getResidueNumber() == info_b->getResidueNumber()
      && info_a->getResidueName() != info_b->getResidueName()) return false;

  const std::string res_a = info_a->getResidueName();
  const std::string res_b = info_b->getResidueName();
  if (!definitions::is_polymer(res_a) || !definitions::is_polymer(res_b)) {
    return false;
  }

  return std::abs(info_a->getResidueNumber() - info_b->getResidueNumber()) <= threshold;
}

ContactRecipe<AtomRec, AtomRec, HBondParams> make_hbond_recipe() {
  return {
    HBondParams{},
    +[](const AtomRec& rec) { return (rec.type & AtomType::HbondDonor)    == AtomType::HbondDonor; },
    +[](const AtomRec& rec) { return (rec.type & AtomType::HbondAcceptor) == AtomType::HbondAcceptor; },
    +[](u32 a, u32 b, float d_sq, const ContactContext& ctx) -> InteractionType {
      const auto& opts = ctx.get_params<HBondParams>();
      const auto &donor    = ctx.topology.atom(a).atom.get();
      const auto &acceptor = ctx.topology.atom(b).atom.get();

      double max_dist = (donor.getAtomicNum() == Element::S || acceptor.getAtomicNum() == Element::S)
                            ? opts.max_sulfur_dist
                            : opts.max_dist;

      if (d_sq < 2.0 || d_sq > max_dist * max_dist) return InteractionType::None;
      if (are_molstar_residueids_close(ctx.molecule(), donor, acceptor, 0))  return InteractionType::None;
      if (!opts.include_water && is_water_hbond(donor, acceptor))            return InteractionType::None;
      if (!are_geometrically_viable(ctx.molecule(), ctx.conformer(), donor, acceptor, opts)) return InteractionType::None;

      return InteractionType::HydrogenBond;
    }
  };
}

ContactRecipe<AtomRec, AtomRec, HBondParams> make_weak_hbond_recipe() {
  return {
    HBondParams{},
    +[](const AtomRec& rec) { return (rec.type & AtomType::WeakHbondDonor) == AtomType::WeakHbondDonor; },
    +[](const AtomRec& rec) { return (rec.type & AtomType::HbondAcceptor) == AtomType::HbondAcceptor; },
    +[](u32 a, u32 b, float d_sq, const ContactContext& ctx) -> InteractionType {
      const auto& params = ctx.get_params<HBondParams>();
      const auto &wdonor   = ctx.topology.atom(a).atom;
      const auto &acceptor = ctx.topology.atom(b).atom;

      // FIX: we now need this additional distance check, because we're not passing down parameters!
      // We could make a WeakHBondParameters struct that inherits from HBondParameters and overrides the max_dist.
      if (d_sq > params.max_dist * params.max_dist) return InteractionType::None;
      if (are_molstar_residueids_close(ctx.molecule(), wdonor, acceptor, 0))  return InteractionType::None;
      if (!params.include_water && is_water_hbond(wdonor, acceptor))          return InteractionType::None;
      if (!are_geometrically_viable(ctx.molecule(), ctx.conformer(), wdonor, acceptor, params)) return InteractionType::None;

      return InteractionType::WeakHydrogenBond;
    }
  };
}

} // namespace lahuta::molstar
