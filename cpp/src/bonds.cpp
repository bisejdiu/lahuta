#include "bonds.hpp"
#include "bonds/lookup.hpp"
#include "bonds/table.hpp"
#include "common.hpp"
#include "convert.hpp"
#include <rdkit/GraphMol/PeriodicTable.h>
#include <vector>

namespace lahuta {

// FIX: see if this fixes the issue with the periodic table
static const RDKit::PeriodicTable *tbl = RDKit::PeriodicTable::getTable();

BondAssignmentResult assign_bonds(RDKit::RWMol &mol, const NSResults &results) {

  std::vector<int> non_predef_atom_indices;
  non_predef_atom_indices.reserve(mol.getNumAtoms());
  std::vector<std::pair<int, int>> bonds;
  std::vector<bool> seen(mol.getNumAtoms(), false);

  std::vector<float> rcov;
  rcov.resize(mol.getNumAtoms());
  for (const auto atom : mol.atoms()) {
    rcov[atom->getIdx()] = tbl->getRcovalent(atom->getAtomicNum());
  }

  for (auto i = 0; i < results.get_pairs().size(); i++) {
    const auto &[index_1, index_2] = results.get_pairs()[i];
    auto dist_sq = results.get_distances()[i];

    AtomInfo atom_1(mol, index_1);
    AtomInfo atom_2(mol, index_2);

    if (!atom_1.info || !atom_2.info) continue;
    if (atom_1.is_hydrogen && atom_2.is_hydrogen) continue;
    if (!is_same_conformer(atom_1, atom_2)) continue;

    PossiblyBonded bonded = getIntraBondOrder(atom_1, atom_2);

    if (bonded.bond_type != RDKit::Bond::BondType::UNSPECIFIED) {

      double pair_thr = get_pair_threshold(atom_1.atom->getAtomicNum(), atom_2.atom->getAtomicNum());

      if (dist_sq <= pair_thr * pair_thr) {
        if (atom_1.is_hydrogen ^ atom_2.is_hydrogen) {
          auto non_h_atom = atom_1.is_hydrogen ? atom_2.atom : atom_1.atom;
          // branchless
          // RDKit::Atom* non_h_atom = a + ((b - a) & -(is_a_h));
          non_h_atom->setNumExplicitHs(non_h_atom->getNumExplicitHs() + 1);
        }
        mol.addBond(atom_1.atom->getIdx(), atom_2.atom->getIdx(), bonded.bond_type);
      }

    } else if (!bonded.atom1_is_predef && !bonded.atom2_is_predef) {

      if (!seen[atom_1.atom->getIdx()]) {
        non_predef_atom_indices.push_back(atom_1.atom->getIdx());
        seen[atom_1.atom->getIdx()] = true;
      }

      if (!seen[atom_2.atom->getIdx()]) {
        non_predef_atom_indices.push_back(atom_2.atom->getIdx());
        seen[atom_2.atom->getIdx()] = true;
      }

      if (is_bonded_obmol(atom_1.atom, atom_2.atom, dist_sq, 0.45, rcov)) {
        bonds.emplace_back(atom_1.atom->getIdx(), atom_2.atom->getIdx());
      }
    } else {
      // if we are here then we have a predef atom and a non-predef atom

      gemmi::El e1 = gemmi::find_element(atom_1.atom->getSymbol().c_str());
      gemmi::El e2 = gemmi::find_element(atom_2.atom->getSymbol().c_str());
      if (gemmi::is_metal(e1) || gemmi::is_metal(e2)) continue;

      double pair_thr = get_pair_threshold(atom_1.atom->getAtomicNum(), atom_2.atom->getAtomicNum());
      if (dist_sq <= pair_thr * pair_thr) {
        if (atom_1.is_hydrogen ^ atom_2.is_hydrogen) {
          auto non_h_atom = atom_1.is_hydrogen ? atom_2.atom : atom_1.atom;
          non_h_atom->setNumExplicitHs(non_h_atom->getNumExplicitHs() + 1);
        }
        mol.addBond(atom_1.atom->getIdx(), atom_2.atom->getIdx(), RDKit::Bond::BondType::SINGLE);
      }
    }
  }

  if (non_predef_atom_indices.empty()) return {};

  std::vector<int> index_mapping;
  index_mapping.resize(mol.getNumAtoms(), -1);
  for (size_t i = 0; i < non_predef_atom_indices.size(); ++i) {
    index_mapping[non_predef_atom_indices[i]] = i;
  }

  auto new_mol = filter_with_conf(mol, non_predef_atom_indices);
  for (const auto &bond : bonds) {
    int aIx = index_mapping[bond.first];
    int bIx = index_mapping[bond.second];

    if (!new_mol.getBondBetweenAtoms(aIx, bIx)) {
      new_mol.addBond(aIx, bIx, RDKit::Bond::BondType::SINGLE);
    }
  }

  return {new_mol, non_predef_atom_indices};
};

} // namespace lahuta
