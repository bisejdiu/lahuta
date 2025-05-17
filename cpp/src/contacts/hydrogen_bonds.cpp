#include "contacts/hydrogen_bonds.hpp"
#include "chemistry/neighbors.hpp"
#include "chemistry/predicates.hpp"
#include "common.hpp"
#include "contacts/hbond_geo_validity.hpp"
#include "contacts/search.hpp"
#include "contacts/utils.hpp"
#include "elements.hpp"
#include "lahuta.hpp"

// clang-format off
namespace lahuta {

// FIX: use the new syntax to check for waters
const std::array<std::string, 11> WaterResidues = {"HOH", "W", "SOL", "TIP3", "SPC", "H2O", "TIP4", "TIP", "DOD", "D3O", "WAT"};

bool is_water(const RDKit::Atom &atom) {
  auto res_info = static_cast<const RDKit::AtomPDBResidueInfo *>(atom.getMonomerInfo());
  if (!res_info) return false;
  return common::contains(WaterResidues, res_info->getResidueName());
}

bool is_water_hbond(const RDKit::Atom &atom_a, const RDKit::Atom &atom_b) {
  return is_water(atom_a) && is_water(atom_b);
}

AtomType add_hydrogen_donor(const RDKit::RWMol &mol, const RDKit::Atom &atom) {
  int total_h = atom.getNumExplicitHs() + atom.getNumCompImplicitHs();

  // include both nitrogen atoms in histidine due to their often ambiguous protonation assignment
  if (chemistry::is_histidine_nitrogen(atom, mol)) return AtomType::HBOND_DONOR;

  // nitrogen, oxygen, or sulfur with hydrogen attached
  int atomic_num = atom.getAtomicNum();
  if (total_h > 0 && (atomic_num == Element::N || atomic_num == Element::O || atomic_num == Element::S)) {
    return AtomType::HBOND_DONOR;
  }

  return AtomType::NONE;
}

AtomType add_hydrogen_acceptor(const RDKit::RWMol &mol, const RDKit::Atom &atom) {
  auto *res_info = static_cast<const RDKit::AtomPDBResidueInfo *>(atom.getMonomerInfo());
  if (!res_info) return AtomType::NONE;

  int formal_charge = atom.getFormalCharge();
  int atomic_num = atom.getAtomicNum();

  // Assume all oxygen atoms are acceptors!
  if (atomic_num == Element::O) return AtomType::HBOND_ACCEPTOR;

  if (atomic_num == Element::N) {
    // include both nitrogen atoms in histidine due to their often ambiguous protonation assignment
    if (chemistry::is_histidine_nitrogen(atom, mol)) return AtomType::HBOND_ACCEPTOR;
    if (formal_charge < 1) {
      // Neutral nitrogen might be an acceptor
      // It must have at least one lone pair not conjugated
      unsigned int total_bonds = get_bond_count(mol, atom) + atom.getNumCompImplicitHs();

      auto hybridization = atom.getHybridization();
      if (   (hybridization == HybridizationType::SP3 && total_bonds < 4)
          || (hybridization == HybridizationType::SP2 && total_bonds < 3)
          || (hybridization == HybridizationType::SP  && total_bonds < 2)) {
        return AtomType::HBOND_ACCEPTOR;
      }
    }
  }

  if (atomic_num == Element::S) {
    std::string res_name = res_info->getResidueName();
    if (res_name == "CYS" || res_name == "MET" || formal_charge == -1) return AtomType::HBOND_ACCEPTOR;
  }

  return AtomType::NONE;
}

AtomType add_weak_hydrogen_donor(const RDKit::RWMol &mol, const RDKit::Atom &atom) {
  int total_h = atom.getNumExplicitHs() + atom.getNumCompImplicitHs();

  if (atom.getAtomicNum() == Element::C && total_h > 0) {
    if (get_bond_count(mol, atom, Element::N) > 0 ||
        get_bond_count(mol, atom, Element::O) > 0 ||
        chemistry::in_aromatic_ring_with_N_or_O(mol, atom)) {
      return AtomType::WEAK_HBOND_DONOR;
    }
  }

  return AtomType::NONE;
}

Contacts find_hydrogen_bonds(const Luni &luni, std::optional<HBondParameters> params) {

  const double distFactor = 1.2;
  const constexpr double MAX_LINE_OF_SIGHT_DISTANCE = 3.0;
  HBondParameters opts = params.value_or(HBondParameters{});

  Contacts contacts(&luni); // FIX: Contacts requires the Luni object (remove?).
  auto &mol = luni.get_molecule();

  const auto donors    = AtomEntityCollection::filter(&luni, AtomType::HBOND_DONOR);
  const auto acceptors = AtomEntityCollection::filter(&luni, AtomType::HBOND_ACCEPTOR);

  auto results =
      EntityNeighborSearch::search(donors, acceptors, std::max(opts.max_dist, opts.max_sulfur_dist));

  std::cout << "results: " << results.size() << std::endl;
  std::unordered_set<std::pair<int, int>, common::PairHash> seen;

  const double distMax = distFactor * MAX_LINE_OF_SIGHT_DISTANCE;
  FastNS grid = FastNS(luni.get_conformer().getPositions());
  auto ok = grid.build(distMax);
  if (!ok) {
    std::cerr << "Failed to build the grid" << std::endl;
    return contacts;
  }

  for (const auto &[pair, dist] : results) {
    auto [donor_index, acceptor_index] = pair;
    auto &donor    = donors.get_data()[donor_index].atoms.front();
    auto &acceptor = acceptors.get_data()[acceptor_index].atoms.front();

    // check the line of sight
    auto donor_com = donors.get_data()[donor_index].center;
    auto acceptor_com = acceptors.get_data()[acceptor_index].center;

    RDGeom::Point3D midpoint = (*donor_com + *acceptor_com) / 2;

    auto ns = grid.search({midpoint});

    // bool line_of_sight_blocked = false;
    // for (const auto &[pair, dist] : ns) {
    //   auto &[_, j] = pair;

    //   const RDKit::Atom *atom = mol.getAtomWithIdx(j);

    //   if (atom->getAtomicNum() == 1) continue;

    //   auto vdw = gemmi::vdw_radius(gemmi::El(static_cast<unsigned char>(atom->getAtomicNum())));
    //   if (vdw * vdw * distFactor * distFactor <= dist) continue;

    //   AtomInfo atom_1(mol, donor->getIdx());
    //   AtomInfo atom_2(mol, acceptor->getIdx());
    //   AtomInfo atom_3(mol, atom->getIdx());
    //   if (!is_same_conformer(atom_1, atom_2) || !is_same_conformer(atom_1, atom_3)) {
    //     continue;
    //   }

    //   AtomEntity v = donors.get_data()[donor_index];
    //   AtomEntity w = acceptors.get_data()[acceptor_index];

    //   if (v.has_atom(atom) || w.has_atom(atom)) continue;

    //   auto atom_pos = luni.get_conformer().getAtomPos(atom->getIdx());
    //   if ((compute_dist_sq(*v.center, atom_pos) < 1.0) || (compute_dist_sq(*w.center, atom_pos) < 1.0)) {
    //     continue;
    //   }

    //   line_of_sight_blocked = true;
    // }

    // if (line_of_sight_blocked) continue;

    // FIX: improve this (it only checks residue numbers)
    if (are_residueids_close(mol, *donor, *acceptor, 0)) continue;

    double max_dist = (donor->getAtomicNum() == Element::S || acceptor->getAtomicNum() == Element::S)
                          ? opts.max_sulfur_dist
                          : opts.max_dist;

    if (dist > max_dist * max_dist) continue;

    if (!opts.include_water && is_water_hbond(*donor, *acceptor)) continue;
    if (!hb_geo::are_geometrically_viable(mol, *donor, *acceptor, opts)) continue;

    if (donor->getIdx() == acceptor->getIdx() || dist < 2.0) continue;
    if (is_duplicate({donor->getIdx(), acceptor->getIdx()}, seen)) continue;

    contacts.add(Contact(
        static_cast<EntityID>(donor->getIdx()),
        static_cast<EntityID>(acceptor->getIdx()),
        dist,
        InteractionType::HydrogenBond));
  }

  return contacts;
}

Contacts find_weak_hydrogen_bonds(const Luni &luni, std::optional<HBondParameters> params) {

  Contacts contacts(&luni);
  const auto &mol = luni.get_molecule();
  HBondParameters opts = params.value_or(HBondParameters{});

  const auto weak_donor_atoms = AtomEntityCollection::filter(&luni, AtomType::WEAK_HBOND_DONOR);
  const auto acceptor_atoms = AtomEntityCollection::filter(&luni, AtomType::HBOND_ACCEPTOR);

  auto nbrs = EntityNeighborSearch::search(weak_donor_atoms, acceptor_atoms, opts.max_dist);

  for (const auto &[pair, dist] : nbrs) {
    auto [wdonor_index, acceptor_index] = pair;
    const auto &wdonor = weak_donor_atoms.get_data()[wdonor_index];
    const auto &acceptor = acceptor_atoms.get_data()[acceptor_index];

    if (are_residueids_close(mol, *wdonor.atom, *acceptor.atom, 1)) continue;
    if (!opts.include_water && is_water_hbond(*wdonor.atom, *acceptor.atom)) continue;
    if (!hb_geo::are_geometrically_viable(mol, *wdonor.atom, *acceptor.atom, opts)) continue;

    contacts.add(Contact(
        static_cast<EntityID>(wdonor.atom->getIdx()),
        static_cast<EntityID>(acceptor.atom->getIdx()),
        dist,
        InteractionType::WeakHydrogenBond));
  }

  return contacts;
}


} // namespace lahuta
