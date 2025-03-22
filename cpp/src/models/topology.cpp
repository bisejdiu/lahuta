#include "models/topology.hpp"
#include "atom_types.hpp"
#include "models/fast_lookup.hpp"
#include "models/rings.hpp"
#include "models/ssbonds.hpp"
#include "models/tables.hpp"
#include <vector>

// clang-format off
namespace lahuta {

void read_and_build_model_topology(RDKit::RWMol &mol, RDKit::Conformer &conf, ModelParserResult &P) {

  size_t num_atoms = P.coords.size();
  size_t num_residues = P.get_sequence().size();

  // construct all atom pointers in advance
  AtomInfoContainer atomContainer(num_atoms);
  std::vector<std::unique_ptr<RDKit::Atom>> AtomPtrs;
  AtomPtrs.reserve(num_atoms);

  int atom_idx = 0;
  auto sequence = P.get_sequence(); // NOTE: converts "SEQ" to {'S','E','Q'}
  for (int residue = 0; residue < sequence.size(); ++residue) {
    const auto &aa = sequence[residue];
    const auto &entry = StandardAminoAcidDataTable[aa[0]];
    for (size_t local_atom_index = 0; local_atom_index < entry.size; ++local_atom_index) {
      const char *atom_name = entry.atoms[local_atom_index];
      const int ih = entry.ih[local_atom_index];
      const int at = entry.at[local_atom_index];

      int atom_number = StandardAminoAcidAtomicNumbers[atom_name[0]];

      std::unique_ptr<RDKit::Atom> atom = std::make_unique<RDKit::Atom>(atom_number);
      atom->setMonomerInfo(atomContainer.createAtomInfo(
                entry.atoms[local_atom_index],
                atom_idx + 1,
                entry.name,
                residue + 1).release()
      );

      atom->setNumCompImplicitHs(ih);
      atom->setCompAtomType(at);

      AtomPtrs.push_back(std::move(atom));
      atom_idx++;
    }
  }

  std::vector<std::vector<int>> aromatic_atom_indices;
  std::vector<std::vector<int>> aromatic_bond_indices;
  std::vector<int> sulphur_atom_indices;

  aromatic_atom_indices.reserve(sequence.size() / 7);
  aromatic_bond_indices.reserve(sequence.size() / 7);
  sulphur_atom_indices .reserve(64);

  int max_arom_atom_idx, max_arom_bond_idx; // help allocation when populating ring infos

  mol.preAllocateAtoms(num_atoms);
  conf.reserve(num_atoms);

  atom_idx = 0;
  for (int residue_idx = 0; residue_idx < sequence.size(); ++residue_idx) {
    const auto &aa = sequence[residue_idx];
    const auto &entry = StandardAminoAcidDataTable[aa[0]];

    for (size_t local_atom_index = 0; local_atom_index < entry.size; ++local_atom_index) {

      const char* atom_name = entry.atoms[local_atom_index];
      int atom_number = StandardAminoAcidAtomicNumbers[atom_name[0]];

      if (atom_number == 16 && aa[0] == 'C') { // CYS (not MET)
          sulphur_atom_indices.push_back(atom_idx);
      }
      mol.addAtomToBatch(AtomPtrs[atom_idx].release());

      ++atom_idx;
    }

    int residue_start_idx = atom_idx - entry.size;

    // inline lambda to construct and record an aromatic ring.
    auto add_aromatic_ring = [&](const auto &indices, bool skip_first = false) {

      static const int max_ring_index = 3;
      auto max_bond_index = aa[0] == 'W' || aa[0] == 'H' ? 4 : 5;

      auto ring = make_aromatic_ring(residue_start_idx, indices);

      std::vector<int> bonds;
      // if trp, we skip the shared bond (the first one)
      for (size_t ri = static_cast<int>(skip_first); ri < ring.size() - 1; ++ri) {
        bonds.push_back(mol.addBond(ring[ri], ring[ri + 1], RDKit::Bond::AROMATIC));
      }

      bonds.push_back(mol.addBond(ring.back(), ring.front(), RDKit::Bond::AROMATIC));
      if (skip_first == 1) ring.erase(ring.begin());

      aromatic_atom_indices.reserve(ring.size());
      aromatic_bond_indices.reserve(bonds.size());
      aromatic_atom_indices.push_back(ring);
      aromatic_bond_indices.push_back(bonds);

      max_arom_atom_idx = ring[3];
      max_arom_bond_idx = bonds[max_bond_index];
    };

    switch (aa[0]) {
      case 'F':
        add_aromatic_ring(phe_arom_indices);
        break;
      case 'Y':
        add_aromatic_ring(tyr_arom_indices);
        break;
      case 'W':
        add_aromatic_ring(trp_arom_indices5);
        add_aromatic_ring(trp_arom_indices6, 1);
        break;
      case 'H':
        add_aromatic_ring(his_arom_indices);
        break;
      default:
        break;
    }
  }

  // add final OXT atom
  auto entry_ = StandardAminoAcidDataTable[sequence.back()[0]];
  std::unique_ptr<RDKit::Atom> oxt_atom = std::make_unique<RDKit::Atom>(8);
  oxt_atom->setMonomerInfo(atomContainer.createAtomInfo(
      "OXT",
      atom_idx + 1,
      entry_.name,
      sequence.size() + 1).release());

  // add OXT, and set the number of implicit Hs for OXT and the first N atom
  oxt_atom->setNumCompImplicitHs(0);
  oxt_atom->setCompAtomType(9217);
  mol.getAtomWithIdx(0)->setNumCompImplicitHs(3);
  mol.addAtomToBatch(oxt_atom.release());

  // set all atom positions
  conf.setAllAtomPositions(std::move(P.coords));
  mol.addConformer(&conf, true);

  // add bonds
  int residue_start_idx = 0;
  for (int residue = 0; residue < sequence.size(); ++residue) {
    const auto &aa = sequence[residue];

    char aa_type = aa[0];
    const auto &entry = StandardAminoAcidDataTable[aa_type];
    const auto &edges = StandardAminoAcidBondTable[aa_type];

    for (size_t j = 0; j < edges.size; ++j) {

      const auto &edge = edges.edges[j];
      if (edge.order == BondType::AROMATIC) continue;

      int atom1_idx = residue_start_idx + edge.i;
      int atom2_idx = residue_start_idx + edge.j;

      mol.addBond(atom1_idx, atom2_idx, edge.order);
    }

    // peptide bond
    if (residue < sequence.size() - 1) {               // not valid for the last residue
      int c_atom_idx = residue_start_idx + 2;          // C atom is at index 2
      int n_atom_idx = residue_start_idx + entry.size; // N atom of next residue
      mol.addBond(c_atom_idx, n_atom_idx, RDKit::Bond::SINGLE);
    } else {
      // For the last residue, add bond to OXT
      int c_atom_idx = residue_start_idx + 2;
      int oxt_atom_idx = atom_idx;                      // OXT is the last atom added
      mol.addBond(c_atom_idx, oxt_atom_idx, RDKit::Bond::SINGLE);
    }

    residue_start_idx += entry.size;
  }

  // handle disulfide bonds, and bonded S atoms
  auto disulfide_pairs = find_disulfide_bonds(sulphur_atom_indices, conf.getPositions());
  for (const auto &pair : disulfide_pairs) {
    mol.addBond(pair.first, pair.second, BondType::SINGLE);

    // correct the implicit Hs from S atoms
    mol.getAtomWithIdx(pair.first )->setNumCompImplicitHs(0);
    mol.getAtomWithIdx(pair.second)->setNumCompImplicitHs(0);

    // correct the atom type for S atoms
    auto first_at  = static_cast<AtomType>(mol.getAtomWithIdx(pair.first )->getCompAtomType());
    auto second_at = static_cast<AtomType>(mol.getAtomWithIdx(pair.second)->getCompAtomType());
    mol.getAtomWithIdx(pair.first )->setCompAtomType(static_cast<int>(first_at  ^ AtomType::HBOND_DONOR));
    mol.getAtomWithIdx(pair.second)->setCompAtomType(static_cast<int>(second_at ^ AtomType::HBOND_DONOR));
  }

  // Finally, we'll handle aromatic rings
  if (mol.getRingInfo()->isInitialized()) {
    mol.getRingInfo()->reset();
  }
  mol.getRingInfo()->initialize(RDKit::FIND_RING_TYPE_SYMM_SSSR);
  mol.getRingInfo()->addAllRings(aromatic_atom_indices, aromatic_bond_indices, max_arom_atom_idx, max_arom_bond_idx);

  // FIX: we may need to assign the underlying data holders to non-sentinel values
  // if skipping this step leads to any issues
  /*mol.updatePropertyCache(false);*/
}


} // namespace lahuta
