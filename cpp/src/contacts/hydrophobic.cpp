#include "contacts/hydrophobic.hpp"
#include "lahuta.hpp"

namespace lahuta {

AtomType add_hydrophobic_atom(const RDKit::RWMol &mol, const RDKit::Atom &atom) {
  const int atomic_num = atom.getAtomicNum();
  bool flag = false;

  if (atomic_num == 6) {
    flag = true;

    for (const auto &bond : mol.atomBonds(&atom)) {
      const RDKit::Atom *neighbor = bond->getOtherAtom(&atom);
      int neighbor_atomic_num = neighbor->getAtomicNum();

      if (neighbor_atomic_num != 6 && neighbor_atomic_num != 1) {
        flag = false;
        break;
      }
    }
  } else if (atomic_num == 9) { // Fluorine atom
    flag = true;
  }

  if (flag) {
    return AtomType::HYDROPHOBIC;
  }
  return AtomType::NONE;
}

void find_hydrophobic_bonds(Luni &luni, const GeometryOptions &opts, Contacts &container) {

  auto max_dist_sq = 4.0 * 4.0;
  const auto &mol = luni.get_molecule();
  const auto hydrophobic_atoms = get_atom_data(&luni, AtomType::HYDROPHOBIC);

  auto grid = FastNS(hydrophobic_atoms.positions(), std::sqrt(max_dist_sq));
  auto nbrs = grid.self_search();

  for (const auto &[pair, dist] : nbrs) {
    auto [atom1_index, atom2_index] = pair;
    const auto &atom1_data = hydrophobic_atoms.data[atom1_index];
    const auto &atom2_data = hydrophobic_atoms.data[atom2_index];

    // FIX: if `self_search` indeed handles self-pairs, this check is redundant
    if (atom1_data.atom->getIdx() == atom2_data.atom->getIdx()) {
      continue;
    }

    if (are_residueids_close(mol, *atom1_data.atom, *atom2_data.atom, 0)) {
      continue;
    }

    // No hydrophobic interaction if both are fluorine
    if (atom1_data.atom->getAtomicNum() == 9 && atom2_data.atom->getAtomicNum() == 9) {
      continue;
    }

    container.add(Contact(
        EntityID(atom1_data.atom->getIdx()),
        EntityID(atom2_data.atom->getIdx()),
        dist,
        InteractionType::Hydrophobic));
  }
}

} // namespace lahuta
