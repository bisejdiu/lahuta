#include "chemistry/neighbors.hpp"
#include "chemistry/predicates.hpp"
#include "chemistry/types/hbonding.hpp"
#include "chemistry/utils.hpp"
#include "chemistry/elements.hpp"

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

} // namespace lahuta
