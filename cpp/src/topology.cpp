#include "topology.hpp"
#include "bond_order.hpp"
#include "bonds.hpp"
#include "ob/clean_mol.hpp"

// clang-format off
namespace lahuta {

RingEntityCollection Topology::populate_ring_entities() {
    RingEntityCollection rings;
    size_t id = 0;
    for (const std::vector<int> &ring : mol_->getRingInfo()->atomRings()) { // NOTE: at this point, we've already added the rings to the molecule
      rings.add_data(*mol_, ring, id++);
    }
    return rings;
  }

void Topology::compute_bonds(const NSResults &neighbors) {
  BondAssignmentResult result = assign_bonds(*mol_, neighbors);

  // FIX: How do we handle connection records? If we define and use topology option, then
  // we can make this user configurable

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

  mol_->updatePropertyCache(false);
  /*clean_bonds(mol, mol.getConformer());*/
  RDKit::MolOps::setHybridization(*mol_);
  cleanup_predef(*mol_);

  if (result.has_unlisted_resnames) {
    clean_bonds(result.mol, result.mol.getConformer());
    perceive_bond_orders_obabel(result.mol);
    cleanup(result.mol);

    merge_bonds(*mol_, result.mol, result.atom_indices);
  }
}

void Topology::cleanup_predef(RWMol &mol) {
  ROMol::VERTEX_ITER at_begin, at_end;
  boost::tie(at_begin, at_end) = mol.getVertices();
  while (at_begin != at_end) {
    RDKit::Atom *atom = mol[*at_begin];
    atom->calcExplicitValence(false);

    // correct four-valent neutral N -> N+
    // This was github #1029
    if (atom->getAtomicNum() == 7 && atom->getFormalCharge() == 0 && atom->getExplicitValence() == 4) {
      atom->setFormalCharge(1);
    }
    ++at_begin;
  }
}

// FIX: this is called for the parts of the molecule that are not predefined
// It needs to be synchronized with the aromatic residue perception module. There we also get "unknown" residues and process them
// it needs to be clear what's the difference and we need to document it.
void Topology::cleanup(RDKit::RWMol &mol) {
  cleanup_predef(mol);

  // FIXME: provide the mechanism as an option?
  // MolOps::fastFindRings(mol);
  // MolOps::findSSSR(mol);
  bool include_dative_bonds = true;
  MolOps::symmetrizeSSSR(mol, include_dative_bonds);
  MolOps::setAromaticity(mol);
}

void Topology::merge_bonds(RDKit::RWMol &target, RDKit::RWMol &source, const std::vector<int> &index_map) {
  for (auto bondIt = source.beginBonds(); bondIt != source.endBonds(); ++bondIt) {
    const RDKit::Bond *bond = *bondIt;
    int bIdx = index_map[bond->getBeginAtomIdx()];
    int eIdx = index_map[bond->getEndAtomIdx()];
    if (target.getBondBetweenAtoms(bIdx, eIdx) == nullptr) {

      // FIX: todo: instead of doing this check, we should iterate once over
      // all atoms and set the number of explicit hydrogens based on the
      // number of explicitly bonded hydrogens
      auto a = target.getAtomWithIdx(bIdx);
      auto b = target.getAtomWithIdx(eIdx);
      int is_a_h = a->getAtomicNum() == 1;
      int is_b_h = b->getAtomicNum() == 1;

      if (is_a_h ^ is_b_h) {
        auto non_h_atom = a->getAtomicNum() == 1 ? b : a;
        non_h_atom->setNumExplicitHs(non_h_atom->getNumExplicitHs() + 1);
      }
      target.addBond(bIdx, eIdx, bond->getBondType());
    }
  }
}

bool Topology::should_initialize_ringinfo(int mol_size) {
    constexpr int small_threshold = 20'000;
    constexpr int medium_threshold = 50'000;
    constexpr int large_threshold = 100'000;

    if (mol_size < small_threshold) return true;

    // FIX: explain what "filtered" means
    if (mol_size < medium_threshold) {
      spdlog::warn("Filtered molecule size ({}) is large. Performance may be affected.", mol_size);
      return true;
    } else if (mol_size < large_threshold) {
      spdlog::warn("Filtered molecule size ({}) is very large. Performance may be severely affected.", mol_size);
      return true;
    }

    spdlog::warn("Filtered molecule size ({}) is too large. Ring perception will be skipped!", mol_size);
    return false;
  }


} // namespace lahuta
