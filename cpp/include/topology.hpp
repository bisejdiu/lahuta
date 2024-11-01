#ifndef LAHUTA_TOPOLOGY_HPP
#define LAHUTA_TOPOLOGY_HPP

#include "atom_types.hpp"
#include "bond_order.hpp"
#include "bonds.hpp"
#include "contacts/charges.hpp"
#include "contacts/contacts.hpp"
#include "convert.hpp"
#include "ob/clean_mol.hpp"
#include "residues.hpp"
#include "rings.hpp"
#include <rdkit/GraphMol/BondIterators.h>

namespace lahuta {

class Topology {
public:
  std::vector<AtomType> atom_types;
  RingDataVec rings_vec;
  const RDKit::RWMol *mol;
  const Residues *residues;

  Topology(const RDKit::RWMol &mol) : mol(&mol) {}
  Topology() = default;
  ~Topology() { delete residues; }

public:
  void build_residues(const RDKit::RWMol &mol) { residues = new Residues(mol); }

  void assign_arpeggio_atom_types() {

    atom_types.resize(mol->getNumAtoms(), AtomType::NONE);

    for (auto atom : mol->atoms()) {
      AtomType atom_type = get_atom_type(atom);
      atom_types[atom->getIdx()] = atom_type;
    }

    auto unk_indices = residue_props::get_unknown_residues<std::vector<int>>(*residues, AminoAcidNames);
    if (!unk_indices.empty()) {
      std::sort(unk_indices.begin(), unk_indices.end());
      auto new_mol = filter_with_bonds(*mol, unk_indices);
      if (should_initialize_ringinfo(new_mol.getNumAtoms())) {
        auto vec = match_atom_types(new_mol);
        for (size_t i = 0; i < unk_indices.size(); ++i) {
          atom_types[unk_indices[i]] = vec[i];
        }
      }
    }

    rings_vec = create_ringdatavec();
  }

  bool should_initialize_ringinfo(int mol_size) const {
    constexpr int small_threshold = 20'000;
    constexpr int medium_threshold = 50'000;
    constexpr int large_threshold = 100'000;

    if (mol_size < small_threshold) {
      return true;
    } else if (mol_size < medium_threshold) {
      std::cerr << "WARNING: Filtered molecule size (" << mol_size
                << ") is large. Performance may be affected." << std::endl;
      return true;
    } else if (mol_size < large_threshold) {
      std::cerr << "WARNING: Filtered molecule size (" << mol_size
                << ") is very large. Performance may be severely affected." << std::endl;
      return true;
    } else {
      std::cerr << "WARNING: Filtered molecule size (" << mol_size
                << ") is too large. Ring perception will be skipped!" << std::endl;
      return false;
    }
  }

  RingDataVec create_ringdatavec() {
    RingDataVec rings;
    size_t id = 0;
    for (const auto &ring : mol->getRingInfo()->atomRings()) {
      RDGeom::Point3D center, norm;
      RingProps::compute_center(mol, ring, center);
      RingProps::compute_normal(mol, ring, center, norm);

      std::vector<const RDKit::Atom *> atoms;
      for (const int &atom_idx : ring) {
        atoms.push_back(mol->getAtomWithIdx(atom_idx));
      }

      rings.rings.emplace_back(center, norm, atoms, id++);
    }
    return rings;
  }

  void assign_molstar_typing() {

    // FIX: to be replaced by the entitytype manager
    std::vector<AtomType> atom_types_ = AtomTypeAnalysis::analyze(*mol);

    atom_types = std::move(atom_types_);
    rings_vec = create_ringdatavec();
  }

  static void compute_bonds(RDKit::RWMol &mol, const NSResults &neighborResults) {
    auto result = assign_bonds(mol, neighborResults);

    // FIX: Refactor!

    // for (gemmi::Connection &conn : st.connections) {
    //   // Iterate over all models
    //   gemmi::Atom *a1 = st.first_model().find_cra(conn.partner1).atom;
    //   gemmi::Atom *a2 = st.first_model().find_cra(conn.partner2).atom;
    //
    //   if (mol.getBondBetweenAtoms(a1->serial-1, a2->serial-1) == nullptr) {
    //     mol.addBond((unsigned int)a1->serial-1, (unsigned int)a2->serial-1,
    //                 RDKit::Bond::BondType::SINGLE);
    //   }
    // }

    mol.updatePropertyCache(false);
    /*clean_bonds(mol, mol.getConformer());*/
    MolOps::setHybridization(mol);
    cleanup_predef(mol);

    if (result.has_unlisted_resnames) {
      clean_bonds(result.mol, result.mol.getConformer());
      perceive_bond_orders_obabel(result.mol);
      cleanup(result.mol);

      merge_bonds(mol, result.mol, result.atom_indices);
    }
  }

private:
  static void cleanup_predef(RWMol &mol) {
    ROMol::VERTEX_ITER atBegin, atEnd;
    boost::tie(atBegin, atEnd) = mol.getVertices();
    while (atBegin != atEnd) {
      RDKit::Atom *atom = mol[*atBegin];
      atom->calcExplicitValence(false);

      // correct four-valent neutral N -> N+
      // This was github #1029
      if (atom->getAtomicNum() == 7 && atom->getFormalCharge() == 0 && atom->getExplicitValence() == 4) {
        atom->setFormalCharge(1);
      }
      ++atBegin;
    }
  }

  static void cleanup(RDKit::RWMol &mol) {
    cleanup_predef(mol);

    // FIXME: provide the mechanism as an option?
    // MolOps::fastFindRings(mol);
    // MolOps::findSSSR(mol);
    bool include_dative_bonds = true;
    MolOps::symmetrizeSSSR(mol, include_dative_bonds);
    MolOps::setAromaticity(mol);
  }

  static void
  merge_bonds(RDKit::RWMol &targetMol, RDKit::RWMol &sourceMol, const std::vector<int> &indexMap) {
    for (auto bondIt = sourceMol.beginBonds(); bondIt != sourceMol.endBonds(); ++bondIt) {
      const RDKit::Bond *bond = *bondIt;
      int bIdx = indexMap[bond->getBeginAtomIdx()];
      int eIdx = indexMap[bond->getEndAtomIdx()];
      if (targetMol.getBondBetweenAtoms(bIdx, eIdx) == nullptr) {

        // FIX: todo: instead of doing this check, we should iterate once over
        // all atoms and set the number of explicit hydrogens based on the
        // number of explicitly bonded hydrogens
        auto a = targetMol.getAtomWithIdx(bIdx);
        auto b = targetMol.getAtomWithIdx(eIdx);
        int is_a_h = a->getAtomicNum() == 1;
        int is_b_h = b->getAtomicNum() == 1;

        if (is_a_h ^ is_b_h) {
          auto non_h_atom = a->getAtomicNum() == 1 ? b : a;
          non_h_atom->setNumExplicitHs(non_h_atom->getNumExplicitHs() + 1);
        }
        targetMol.addBond(bIdx, eIdx, bond->getBondType());
      }
    }
  }
};

} // namespace lahuta

#endif // LAHUTA_TOPOLOGY_HPP
