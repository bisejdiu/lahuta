#include "models/topology.hpp"
#include "atom_types.hpp"
#include "logging.hpp"
#include "models/factory.hpp"
#include "models/fast_lookup.hpp"
#include "models/pools.hpp"
#include "models/rings.hpp"
#include "models/ssbonds.hpp"
#include "models/tables.hpp"
#include <vector>

// clang-format off
namespace lahuta {

void build_model_topology_def(std::shared_ptr<RDKit::RWMol> &mol, RDKit::Conformer &conf, ModelParserResult &P) {

  mol->skipAtomCleanupDespiteOwnership(true);
  mol->skipBondCleanupDespiteOwnership(true);

  size_t num_atoms = P.coords.size();

  auto* info_pool = PoolFactory<InfoPool>::getFreshPoolForCurrentThread();
  auto* atom_pool = PoolFactory<AtomPool>::getFreshPoolForCurrentThread();
  auto *bond_pool = PoolFactory<BondPool>::getFreshPoolForCurrentThread();

  std::vector<RDKit::Atom*> AtomPtrs;
  AtomPtrs.reserve(num_atoms);

  auto sequence = P.get_sequence(); // NOTE: converts "SEQ" to {'S','E','Q'}

  std::vector<std::vector<int>> aromatic_atom_indices;
  std::vector<std::vector<int>> aromatic_bond_indices;
  std::vector<int> sulphur_atom_indices;

  aromatic_atom_indices.reserve(sequence.size() / 7);
  aromatic_bond_indices.reserve(sequence.size() / 7);
  sulphur_atom_indices .reserve(64);

  std::uint32_t atom_idx = 0;
  for (int residue_idx = 0; residue_idx < sequence.size(); ++residue_idx) {
    const auto &aa = sequence[residue_idx];
    const auto &entry = StandardAminoAcidDataTable[aa[0]];

    for (size_t local_atom_index = 0; local_atom_index < entry.size; ++local_atom_index) {

      const char* atom_name = entry.atoms[local_atom_index];
      int atom_number = StandardAminoAcidAtomicNumbers[atom_name[0]];
      int ih = entry.ih[local_atom_index];
      int at = entry.at[local_atom_index];

      if (atom_number == 16 && aa[0] == 'C') { // CYS (not MET)
          sulphur_atom_indices.push_back(atom_idx);
      }

      RDKit::Atom* atom = atom_pool->createAtom(atom_number);
      atom->setMonomerInfo(info_pool->createAtomInfo(
                entry.atoms[local_atom_index],
                atom_idx + 1,
                entry.name,
                residue_idx + 1)
      );

      atom->setNumCompImplicitHs(ih);
      atom->setCompAtomType(at);

      mol->addAtom(atom, false, true);

      ++atom_idx;
    }

    int residue_start_idx = atom_idx - entry.size;

    // inline lambda to construct and record an aromatic ring.
    auto add_aromatic_ring = [&](const auto &indices, bool skip_first = false) {

      auto ring = make_aromatic_ring(residue_start_idx, indices);

      std::vector<int> bonds;
      // if trp, we skip the shared bond (the first one)
      for (size_t ri = static_cast<int>(skip_first); ri < ring.size() - 1; ++ri) {
        auto bond = bond_pool->createBond(ring[ri], ring[ri + 1], RDKit::Bond::AROMATIC);
        bonds.push_back(mol->addBond(bond, true));
      }

      auto *bond = bond_pool->createBond(ring.back(), ring.front(), RDKit::Bond::AROMATIC);
      bonds.push_back(mol->addBond(bond, true));
      if (skip_first == 1) ring.erase(ring.begin());

      aromatic_atom_indices.reserve(ring.size());
      aromatic_bond_indices.reserve(bonds.size());
      aromatic_atom_indices.push_back(ring);
      aromatic_bond_indices.push_back(bonds);
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
  RDKit::Atom* oxt_atom = atom_pool->createAtom(8);
  oxt_atom->setMonomerInfo(info_pool->createAtomInfo(
      "OXT",
      atom_idx + 1,
      entry_.name,
      sequence.size() + 1));

  // add OXT, and set the number of implicit Hs for OXT and the first N atom
  oxt_atom->setNumCompImplicitHs(0);
  oxt_atom->setCompAtomType(9217);
  mol->getAtomWithIdx(0)->setNumCompImplicitHs(3);
  mol->addAtom(oxt_atom, false, true);

  // set all atom positions
  conf.reserve(num_atoms);
  conf.setAllAtomPositions(std::move(P.coords));
  mol->addConformer(&conf, true);

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

      auto bond = bond_pool->createBond(atom1_idx, atom2_idx, edge.order);
      mol->addBond(bond, true);
    }

    residue_start_idx += entry.size;
  }

  // handle disulfide bonds, and bonded S atoms
  auto disulfide_pairs = find_disulfide_bonds(sulphur_atom_indices, conf.getPositions());
  for (const auto &pair : disulfide_pairs) {
    auto *bond = bond_pool->createBond(pair.first, pair.second, RDKit::Bond::SINGLE);
    mol->addBond(bond, true);

    // correct the implicit Hs from S atoms
    mol->getAtomWithIdx(pair.first )->setNumCompImplicitHs(0);
    mol->getAtomWithIdx(pair.second)->setNumCompImplicitHs(0);

    // correct the atom type for S atoms
    auto first_at  = static_cast<AtomType>(mol->getAtomWithIdx(pair.first )->getCompAtomType());
    auto second_at = static_cast<AtomType>(mol->getAtomWithIdx(pair.second)->getCompAtomType());
    mol->getAtomWithIdx(pair.first )->setCompAtomType(static_cast<int>(first_at  ^ AtomType::HBOND_DONOR));
    mol->getAtomWithIdx(pair.second)->setCompAtomType(static_cast<int>(second_at ^ AtomType::HBOND_DONOR));
  }
  Logger::get_logger()->info("Molecule Atoms: {}", mol->getNumAtoms());
  Logger::get_logger()->info("Molecule Bonds: {}", mol->getNumBonds());

  if (mol->getRingInfo()->isInitialized()) {
    mol->getRingInfo()->reset();
  }

  mol->getRingInfo()->initialize(RDKit::FIND_RING_TYPE_SYMM_SSSR);
  mol->getRingInfo()->addAllRings(aromatic_atom_indices, aromatic_bond_indices);
}

void build_model_topology_csr(std::shared_ptr<RDKit::RWMol> &mol, RDKit::Conformer &conf, ModelParserResult &P) {

  const auto    sequence = P.get_sequence();
  const size_t num_atoms = P.coords.size();
  const size_t num_residues = sequence.size();

  // memory pools
  auto* info_pool = PoolFactory<InfoPool>::getFreshPoolForCurrentThread();
  auto* atom_pool = PoolFactory<AtomPool>::getFreshPoolForCurrentThread();
  auto *bond_pool = PoolFactory<BondPool>::getFreshPoolForCurrentThread();

  // total bond count to pre-allocate memory
  size_t expected_bond_count = 0;
  for (int residue_idx = 0; residue_idx < num_residues; ++residue_idx) {
    const auto &aa = sequence[residue_idx];
    const AminoAcidEdges& edges = StandardAminoAcidBondTable[aa[0]];
    expected_bond_count += edges.size;
  }

  expected_bond_count += num_residues / 10;

  std::vector<RDKit::Atom*> vertices;
  vertices.reserve(num_atoms + 1);  // +1 for OXT

  std::vector<std::vector<int>> aromatic_atom_indices;
  std::vector<std::vector<int>> aromatic_bond_indices;
  std::vector<int> sulphur_atom_indices;

  aromatic_atom_indices.reserve(num_residues / 7);
  aromatic_bond_indices.reserve(num_residues / 7);
  sulphur_atom_indices.reserve(64);

  std::vector<std::pair<size_t, size_t>> edge_list;
  std::vector<RDKit::Bond*> edge_props;
  edge_list.reserve(expected_bond_count);
  edge_props.reserve(expected_bond_count);

  int atom_idx = 0;
  for (int residue_idx = 0; residue_idx < num_residues; ++residue_idx) {
    const auto &aa = sequence[residue_idx];
    const auto &entry = StandardAminoAcidDataTable[aa[0]];

    int residue_start_idx = atom_idx; // residue start index

    // now we create atoms for this residue
    for (size_t local_atom_index = 0; local_atom_index < entry.size; ++local_atom_index) {
      const char *atom_name = entry.atoms[local_atom_index];
      const int ih = entry.ih[local_atom_index];
      const int at = entry.at[local_atom_index];
      int atom_number = StandardAminoAcidAtomicNumbers[atom_name[0]];

      RDKit::Atom* atom = atom_pool->createAtom(atom_number);
      atom->setMonomerInfo(info_pool->createAtomInfo(
                atom_name,
                atom_idx + 1,
                entry.name,
                residue_idx + 1)
      );

      atom->setNumCompImplicitHs(ih);
      atom->setCompAtomType(at);
      atom->setIdx(atom_idx);

      if (atom_number == 16 && aa[0] == 'C') { // CYS sulfur
        sulphur_atom_indices.push_back(atom_idx);
      }

      vertices.push_back(atom);
      atom_idx++;
    }

    const AminoAcidEdges& edges = StandardAminoAcidBondTable[aa[0]];

    int first_bond_idx = edge_props.size();
    for (size_t j = 0; j < edges.size; ++j) {
      const auto& edge = edges.edges[j];
      int atom1_idx = residue_start_idx + edge.i;
      int atom2_idx = residue_start_idx + edge.j;

      auto bond = bond_pool->createBond(atom1_idx, atom2_idx, edge.order);
      bond->setIdx(edge_props.size());

      edge_list.emplace_back(atom1_idx, atom2_idx);
      edge_props.push_back(bond);
    }

    // generic lambda
    auto process_ring = [&aromatic_atom_indices, &aromatic_bond_indices, residue_start_idx, first_bond_idx]
                      (const auto& atom_indices, const auto& bond_indices) {
        std::vector<int> ring_atom_indices;
        std::vector<int> ring_bond_indices;

        // map predefined aromatic atom indices to actual atom indices
        for (int idx : atom_indices) {
            ring_atom_indices.push_back(residue_start_idx + idx);
        }

        // find bond indices
        for (size_t i = 0; i < atom_indices.size(); ++i) {
            ring_bond_indices.push_back(first_bond_idx + bond_indices[i]);
        }

        aromatic_atom_indices.push_back(ring_atom_indices);
        aromatic_bond_indices.push_back(ring_bond_indices);
    };

    switch (aa[0]) {
        case 'F':
            process_ring(phe_arom_indices, phe_bond_indices);
            break;
        case 'Y':
            process_ring(tyr_arom_indices, tyr_bond_indices);
            break;
        case 'H':
            process_ring(his_arom_indices, his_bond_indices);
            break;
        case 'W':
            // Tryptophan has two rings to process
            process_ring(trp_arom_indices5, trp_bond_indices5);
            process_ring(trp_arom_indices6, trp_bond_indices6);
            break;
    }
  }

  // Add terminal OXT atom
  auto last_entry = StandardAminoAcidDataTable[sequence.back()[0]];
  RDKit::Atom* oxt_atom = atom_pool->createAtom(8);
  oxt_atom->setMonomerInfo(info_pool->createAtomInfo(
      "OXT",
      atom_idx + 1,
      last_entry.name,
      num_residues));
  oxt_atom->setNumCompImplicitHs(0);
  oxt_atom->setCompAtomType(9217);
  oxt_atom->setIdx(atom_idx);

  vertices.push_back(oxt_atom);
  // for the first atom, set implicit Hs to 3
  vertices.front()->setNumCompImplicitHs(3);

  // positions
  conf.reserve(num_atoms);
  conf.setAllAtomPositions(std::move(P.coords));

  // disulfide bonds
  auto disulfide_pairs = find_disulfide_bonds(sulphur_atom_indices, conf.getPositions());
  for (const auto &pair : disulfide_pairs) {
    RDKit::Bond* bond = bond_pool->createBond(pair.first, pair.second, RDKit::Bond::SINGLE);
    bond->setIdx(edge_props.size());

    auto it = std::lower_bound(edge_list.begin(), edge_list.end(), std::make_pair(pair.first, pair.second));
    edge_list.insert(it, std::make_pair(pair.first, pair.second));
    edge_props.push_back(bond);
  }

  mol = std::make_shared<RDKit::RWMol>(vertices, edge_list, edge_props, GraphType::CSRMolGraph);
  mol->addConformer(&conf, true);

  for (const auto &pair : disulfide_pairs) {
    // Correct the implicit Hs from S atoms
    mol->getAtomWithIdx(pair.first) ->setNumCompImplicitHs(0);
    mol->getAtomWithIdx(pair.second)->setNumCompImplicitHs(0);

    // Correct the atom type for S-S bound atoms
    auto first_at   = static_cast<AtomType>(mol->getAtomWithIdx(pair.first) ->getCompAtomType());
    auto second_at  = static_cast<AtomType>(mol->getAtomWithIdx(pair.second)->getCompAtomType());
    mol->getAtomWithIdx(pair.first) ->setCompAtomType(static_cast<int>(first_at  ^ AtomType::HBOND_DONOR));
    mol->getAtomWithIdx(pair.second)->setCompAtomType(static_cast<int>(second_at ^ AtomType::HBOND_DONOR));
  }

  if (mol->getRingInfo()->isInitialized()) {
    mol->getRingInfo()->reset();
  }
  mol->getRingInfo()->initialize(RDKit::FIND_RING_TYPE_SYMM_SSSR);
  mol->getRingInfo()->addAllRings(aromatic_atom_indices, aromatic_bond_indices);
}

} // namespace lahuta
