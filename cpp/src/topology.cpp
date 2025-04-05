#include "topology.hpp"
#include "aromatics.hpp"
#include "bond_order.hpp"
#include "bonds.hpp"
#include "convert.hpp"
#include "definitions.hpp"
#include "logging.hpp"
#include "ob/clean_mol.hpp"

// clang-format off
namespace lahuta {

void Topology::build(TopologyBuildingOptions tops) {

  if (!mol_) {
    Logger::get_logger()->critical("Cannot build topology without a molecule.");
    throw std::runtime_error("Make sure to provide a molecule before building the topology.");
  }

  try {
    if (tops.compute_bonds) {
      auto start = std::chrono::high_resolution_clock::now();
      auto grid = FastNS(mol_->getConformer().getPositions());
      auto ok = grid.build(tops.cutoff);
      if (!ok) {
        throw std::runtime_error("Failed to build the grid for neighbor search.");
      }
      auto neighbors = std::make_shared<NSResults>(grid.self_search());

      // FIX: neighbor computation can technically be the responsibility of compute_bonds, but
      // that moves them perhaps too much down the stack, and makes control of the process
      // more difficult (e.g., if we need to use a memory pool or arena allocator)
      this->compute_bonds(*neighbors);
      auto end = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
      Logger::get_logger()->info("Computed bonds in {} ms", duration.count());
    }

    // build residue information
    residues->build();

    // populate ring information to RDKit Mol
    auto start = std::chrono::high_resolution_clock::now();
    initialize_and_populate_ringinfo(*mol_, *residues);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    Logger::get_logger()->info("Populated ring information in {} us", duration.count());

    auto start2 = std::chrono::high_resolution_clock::now();
    switch (tops.atom_typing_method) {
      case ContactComputerType::Molstar:
        this->assign_molstar_typing();
        break;
      case ContactComputerType::Arpeggio:
        this->assign_arpeggio_atom_types();
        break;
      default:
        break;
    }
    auto end2 = std::chrono::high_resolution_clock::now();
    auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(end2 - start2);
    Logger::get_logger()->info("Assigned molstar types in {} us", duration2.count());

  } catch (const std::runtime_error &e) {
    Logger::get_logger()->critical(
        "Error creating topology! Exception caught: {}. Will not terminate, "
        "but no topology-based features will be available.",
        e.what());
  }
}

void Topology::assign_arpeggio_atom_types() {

  atom_types.reserve(mol_->getNumAtoms());
  for (auto atom : mol_->atoms()) {
    AtomType atom_type = get_atom_type(atom);
    atom_types.add_data(*mol_, atom, atom_type);
  }

  auto unk_indices = residues->filter(std::not_fn(definitions::is_protein_extended)).get_atom_ids();
  if (!unk_indices.empty()) {
    std::sort(unk_indices.begin(), unk_indices.end());
    auto new_mol = filter_with_bonds(*mol_, unk_indices);
    if (should_initialize_ringinfo(new_mol.getNumAtoms())) {
      auto vec = match_atom_types(new_mol);
      for (size_t i = 0; i < unk_indices.size(); ++i) {
        atom_types.add_data(*mol_, mol_->getAtomWithIdx(unk_indices[i]), vec[i]);
      }
    }
  }

  rings_vec = populate_ring_entities();
}

RingEntityCollection Topology::populate_ring_entities() {
  RingEntityCollection rings;
  size_t id = 0;
  for (const std::vector<int> &ring :
       // NOTE: at this point, we should have already added the rings to the molecule
       mol_->getRingInfo()->atomRings()) {
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
  for (auto atom : mol.atoms()) {
    atom->calcExplicitValence(false);
    // correct four-valent neutral N -> N+
    // This was github #1029
    if (atom->getAtomicNum() == 7 &&
        atom->getFormalCharge() == 0 &&
        atom->getExplicitValence() == 4) {
      atom->setFormalCharge(1);
    }
  }
}


// FIX: this is called for the parts of the molecule that are not predefined
// It needs to be synchronized with the aromatic residue perception module. There we also get "unknown"
// residues and process them it needs to be clear what's the difference and we need to document it.
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
  for (const auto &bond : source.bonds()) {
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

size_t Topology::total_size() const {
  size_t total = sizeof(*this);

  total += sizeof(AtomEntity) * atom_types.get_data().size();
  total += sizeof(RingEntity) * rings_vec.get_data().size();
  total += sizeof(GroupEntity) * features.get_data().size();

  if (mol_) {
    total += sizeof(*mol_);
  }

  total += residues->total_size();

  return total;
}

bool Topology::should_initialize_ringinfo(int mol_size) {
  constexpr int small_threshold = 20'000;
  constexpr int medium_threshold = 50'000;
  constexpr int large_threshold = 100'000;

  if (mol_size < small_threshold) return true;

  // FIX: explain what "filtered" means
  if (mol_size < medium_threshold) {
    Logger::get_logger()->warn(
        "Filtered molecule size ({}) is large. Performance may be affected.",
        mol_size);
    return true;
  } else if (mol_size < large_threshold) {
    Logger::get_logger()->warn(
        "Filtered molecule size ({}) is very large. Performance may be severely affected.",
        mol_size);
    return true;
  }

  Logger::get_logger()->error(
      "Filtered molecule size ({}) is too large. Ring perception will be skipped!",
      mol_size);
  return false;
}

} // namespace lahuta
