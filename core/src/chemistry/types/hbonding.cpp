#include "chemistry/types/hbonding.hpp"
#include "chemistry/neighbors.hpp"
#include "chemistry/predicates.hpp"
#include "chemistry/utils.hpp"
#include "elements.hpp"

// clang-format off
namespace lahuta {

AtomType add_hydrogen_donor(const RDKit::RWMol &mol, const RDKit::Atom &atom) {
  int total_h = atom.getNumExplicitHs() + atom.getNumCompImplicitHs();

  // include both nitrogen atoms in histidine due to their often ambiguous protonation assignment
  if (chemistry::is_histidine_nitrogen(atom, mol)) return AtomType::HbondDonor;

  // nitrogen, oxygen, or sulfur with hydrogen attached
  int atomic_num = atom.getAtomicNum();
  if (total_h > 0 && (atomic_num == Element::N || atomic_num == Element::O || atomic_num == Element::S)) {
    return AtomType::HbondDonor;
  }

  return AtomType::None;
}

AtomType add_hydrogen_acceptor(const RDKit::RWMol &mol, const RDKit::Atom &atom) {
  auto *res_info = static_cast<const RDKit::AtomPDBResidueInfo *>(atom.getMonomerInfo());
  if (!res_info) return AtomType::None;

  int formal_charge = atom.getFormalCharge();
  int atomic_num = atom.getAtomicNum();

  // Assume all oxygen atoms are acceptors!
  if (atomic_num == Element::O) return AtomType::HbondAcceptor;

  if (atomic_num == Element::N) {
    // include both nitrogen atoms in histidine due to their often ambiguous protonation assignment
    if (chemistry::is_histidine_nitrogen(atom, mol)) return AtomType::HbondAcceptor;
    if (formal_charge < 1) {
      // Neutral nitrogen might be an acceptor
      // It must have at least one lone pair not conjugated
      unsigned int total_bonds = get_bond_count(mol, atom) + atom.getNumCompImplicitHs();

      auto hybridization = atom.getHybridization();
      if (   (hybridization == HybridizationType::SP3 && total_bonds < 4)
          || (hybridization == HybridizationType::SP2 && total_bonds < 3)
          || (hybridization == HybridizationType::SP  && total_bonds < 2)) {
        return AtomType::HbondAcceptor;
      }
    }
  }

  if (atomic_num == Element::S) {
    std::string res_name = res_info->getResidueName();
    if (res_name == "CYS" || res_name == "MET" || formal_charge == -1) return AtomType::HbondAcceptor;
  }

  return AtomType::None;
}

AtomType add_weak_hydrogen_donor(const RDKit::RWMol &mol, const RDKit::Atom &atom) {
  int total_h = atom.getNumExplicitHs() + atom.getNumCompImplicitHs();

  if (atom.getAtomicNum() == Element::C && total_h > 0) {
    if (get_bond_count(mol, atom, Element::N) > 0 ||
        get_bond_count(mol, atom, Element::O) > 0 ||
        chemistry::in_aromatic_ring_with_N_or_O(mol, atom)) {
      return AtomType::WeakHbondDonor;
    }
  }

  return AtomType::None;
}

// ContactSet find_hydrogen_bonds(const Topology &topology, const HBondParameters &opts) {
//   return find_contacts(
//     {topology, opts},
//     [](const AtomRec& rec) { return (rec.type & AtomType::HbondDonor)    == AtomType::HbondDonor; },
//     [](const AtomRec& rec) { return (rec.type & AtomType::HbondAcceptor) == AtomType::HbondAcceptor; },
//     {std::max(opts.max_dist, opts.max_sulfur_dist), 0.7},
//     [](std::uint32_t rec_idx_a, std::uint32_t rec_idx_b, float dist, const ContactContext& ctx) -> InteractionType {
//       const auto& opts = ctx.get_params<HBondParameters>();
//       const auto &donor    = ctx.topology.atom(rec_idx_a).atom;
//       const auto &acceptor = ctx.topology.atom(rec_idx_b).atom;
// 
//       double max_dist = (donor.getAtomicNum() == Element::S || acceptor.getAtomicNum() == Element::S)
//                             ? opts.max_sulfur_dist
//                             : opts.max_dist;
// 
//       if (dist < 2.0 || dist > max_dist * max_dist)        return InteractionType::None;
//       if (are_residueids_close(ctx.molecule(), donor, acceptor, 0)) return InteractionType::None;
//       if (!opts.include_water && is_water_hbond(donor, acceptor))        return InteractionType::None;
//       if (!hb_geo::are_geometrically_viable(ctx.molecule(), donor, acceptor, opts)) return InteractionType::None;
// 
//       return InteractionType::HydrogenBond;
//     }
//   );
// }
// 
// ContactSet find_weak_hydrogen_bonds(const Topology &topology, HBondParameters params) {
//   return find_contacts(
//     {topology, params},
//     [](const AtomRec& rec) { return (rec.type & AtomType::WeakHbondDonor) == AtomType::WeakHbondDonor; },
//     [](const AtomRec& rec) { return (rec.type & AtomType::HbondAcceptor)   == AtomType::HbondAcceptor; },
//     {params.max_dist, 0.10},
//     [](std::uint32_t rec_idx_a, std::uint32_t rec_idx_b, float dist, const ContactContext& ctx) -> InteractionType {
//       const auto& params = ctx.get_params<HBondParameters>();
//       const auto &wdonor   = ctx.topology.atom(rec_idx_a).atom;
//       const auto &acceptor = ctx.topology.atom(rec_idx_b).atom;
// 
//       if (are_residueids_close(ctx.molecule(), wdonor, acceptor, 1))                   return InteractionType::None;
//       if (!params.include_water && is_water_hbond(wdonor, acceptor))        return InteractionType::None;
//       if (!hb_geo::are_geometrically_viable(ctx.molecule(), wdonor, acceptor, params)) return InteractionType::None;
// 
//       return InteractionType::WeakHydrogenBond;
//     }
//   );
// }

} // namespace lahuta
