#include "models/topology.hpp"
#include "typing/types.hpp"
#include "models/factory.hpp"
#include "models/fast_lookup.hpp"
#include "models/pools.hpp"
#include "models/ssbonds.hpp"
#include "models/tables.hpp"
#include <logging.hpp>
#include <vector>

// clang-format off
namespace lahuta {

void build_model_topology_def(std::shared_ptr<RDKit::RWMol> &mol, RDKit::Conformer &conf, const ModelParserResult &P) {

  const auto    sequence     = P.get_sequence();
  const size_t  num_atoms    = P.coords.size();
  const size_t  num_residues = sequence.size();

  auto *info_pool = PoolFactory<InfoPool>::getFreshPoolForCurrentThread();
  auto *atom_pool = PoolFactory<AtomPool>::getFreshPoolForCurrentThread();
  auto *bond_pool = PoolFactory<BondPool>::getFreshPoolForCurrentThread();

  size_t expected_bond_count = 0;
  for (int residue_idx = 0; residue_idx < static_cast<int>(num_residues); ++residue_idx) {
    const auto &aa = sequence[residue_idx];
    const AminoAcidEdges &edges = StandardAminoAcidBondTable[aa[0]];
    expected_bond_count += edges.size + 1; // +1 for peptide bond to next residue
  }

  std::vector<RDKit::Atom*> vertices;
  vertices.reserve(num_atoms); // includes OXT

  std::vector<std::vector<int>> aromatic_atom_indices;
  std::vector<std::vector<int>> aromatic_bond_indices;
  std::vector<int> sulphur_atom_indices;

  std::vector<RDKit::Bond*> bonds;
  bonds.reserve(expected_bond_count);

  int atom_idx = 0;
  for (int residue_idx = 0; residue_idx < static_cast<int>(num_residues); ++residue_idx) {
    const auto &aa    = sequence[residue_idx];
    const auto &entry = StandardAminoAcidDataTable[aa[0]];

    const int residue_start_idx = atom_idx;

    // Create atoms for this residue
    for (size_t local_atom_index = 0; local_atom_index < entry.size; ++local_atom_index) {
      const char *atom_name = entry.atoms[local_atom_index];
      const int   ih        = entry.ih[local_atom_index];
      const int   at        = entry.at[local_atom_index];
      const int   hyb       = entry.hyb[local_atom_index];
      const int   atom_num  = StandardAminoAcidAtomicNumbers[atom_name[0]];

      auto *atom = atom_pool->createAtom();
      atom->setIdx(static_cast<unsigned int>(atom_idx));
      atom->setAtomicNum(atom_num);
      atom->setMonomerInfo(info_pool->createAtomInfo(
          atom_name,
          atom_idx + 1,
          entry.name,
          residue_idx + 1));
      atom->setNumCompImplicitHs(ih);
      atom->setCompAtomType(at);
      atom->setHybridization(static_cast<RDKit::Atom::HybridizationType>(hyb));
      atom->setIdx(static_cast<unsigned int>(atom_idx));

      if (atom_num == 16 && aa[0] == 'C') { // sulfur in CYS only
        sulphur_atom_indices.push_back(atom_idx);
      }

      vertices.push_back(atom);
      ++atom_idx;
    }

    // Residue-internal bonds
    const AminoAcidEdges &edges = StandardAminoAcidBondTable[aa[0]];
    for (size_t j = 0; j < edges.size; ++j) {
      const auto &edge = edges.edges[j];
      const int atom1_idx_local = residue_start_idx + edge.i;
      const int atom2_idx_local = residue_start_idx + edge.j;
      auto bond = bond_pool->createBond(atom1_idx_local, atom2_idx_local, edge.order);
      bonds.push_back(bond);
    }

    // Peptide bond: current C (pos 2) to next residue N (pos 0)
    if (residue_idx < static_cast<int>(num_residues) - 1) {
      const int current_c_idx = residue_start_idx + 2;
      const int next_n_idx    = atom_idx; // start index of next residue
      auto bond = bond_pool->createBond(current_c_idx, next_n_idx, RDKit::Bond::SINGLE);
      bonds.push_back(bond);
    }

    auto process_ring = [&aromatic_atom_indices, &vertices, residue_start_idx](const auto &atom_indices) {
      std::vector<int> ring_atom_indices;
      ring_atom_indices.reserve(atom_indices.size());
      for (int idx : atom_indices) {
        const int global_idx = residue_start_idx + idx;
        ring_atom_indices.push_back(global_idx);
        vertices[static_cast<size_t>(global_idx)]->setIsAromatic(true);
      }
      aromatic_atom_indices.push_back(std::move(ring_atom_indices));
    };

    switch (aa[0]) {
      case 'F': process_ring(phe_arom_indices); break;
      case 'Y': process_ring(tyr_arom_indices); break;
      case 'H': process_ring(his_arom_indices); break;
      case 'W':
          process_ring(trp_arom_indices5);
          process_ring(trp_arom_indices6);
        break;
    }
  }

  // Terminal OXT atom
  if (!sequence.empty()) {
    const auto last_entry = StandardAminoAcidDataTable[sequence.back()[0]];
    RDKit::Atom *oxt_atom = atom_pool->createAtom();
    oxt_atom->setAtomicNum(8);
    oxt_atom->setIdx(static_cast<unsigned int>(atom_idx));
    oxt_atom->setMonomerInfo(info_pool->createAtomInfo(
        "OXT",
        atom_idx + 1,
        last_entry.name,
        static_cast<int>(num_residues)));
    oxt_atom->setNumCompImplicitHs(0);
    oxt_atom->setCompAtomType(9217);
    oxt_atom->setHybridization(RDKit::Atom::SP2);
    oxt_atom->setIdx(static_cast<unsigned int>(atom_idx));
    vertices.push_back(oxt_atom);

    // First atom (N of first residue): tetrahedral geometry and 3 implicit Hs
    if (!vertices.empty()) {
      vertices.front()->setHybridization(RDKit::Atom::SP3);
      vertices.front()->setNumCompImplicitHs(3);
    }

    // C-OXT bond (carbonyl carbon of last residue to OXT)
    const int last_residue_c_idx = atom_idx - static_cast<int>(last_entry.size) + 2;
    auto c_oxt_bond = bond_pool->createBond(last_residue_c_idx, atom_idx, RDKit::Bond::SINGLE);
    bonds.push_back(c_oxt_bond);
  }

  // positions
  conf.reserve(num_atoms);
  conf.setAllAtomPositions(std::move(P.coords));

  // disulfide bonds to edge list
  const auto disulfide_pairs = find_disulfide_bonds(sulphur_atom_indices, conf.getPositions());
  for (const auto &pair : disulfide_pairs) {
    auto bond = bond_pool->createBond(pair.first, pair.second, RDKit::Bond::SINGLE);
    bonds.push_back(bond);
  }

  // Assign stable indices to bonds prior to build
  for (size_t i = 0; i < bonds.size(); ++i) {
    if (bonds[i]) bonds[i]->setIdx(static_cast<unsigned int>(i));
  }

  // build molecule
  mol = std::make_shared<RDKit::RWMol>(vertices, bonds, GraphType::CSRMolGraph);
  mol->addConformer(&conf, true);

  // Build ring bond indices based on final graph
  aromatic_bond_indices.clear();
  aromatic_bond_indices.reserve(aromatic_atom_indices.size());
  for (const auto &ring_atoms : aromatic_atom_indices) {
    std::vector<int> ring_bonds;
    const size_t n = ring_atoms.size();
    if (n < 3) { aromatic_bond_indices.push_back(ring_bonds); continue; }
    for (size_t i = 0; i < n; ++i) {
      const int u = ring_atoms[i];
      const int v = ring_atoms[(i + 1) % n];
      const RDKit::Bond *b = mol->getBondBetweenAtoms(u, v);
      if (b) {
        ring_bonds.push_back(static_cast<int>(b->getIdx()));
      } else {
        Logger::get_logger()->warn("build_model_topology_def: No bond found between atoms {} and {} in ring", u, v);
        ring_bonds.push_back(-1);
      }
    }
    aromatic_bond_indices.push_back(std::move(ring_bonds));
  }

  // Correct properties for S-S bonded atoms now that the graph exists
  for (const auto &pair : disulfide_pairs) {
    mol->getAtomWithIdx(pair.first )->setNumCompImplicitHs(0);
    mol->getAtomWithIdx(pair.second)->setNumCompImplicitHs(0);

    auto first_at  = static_cast<AtomType>(mol->getAtomWithIdx(pair.first )->getCompAtomType());
    auto second_at = static_cast<AtomType>(mol->getAtomWithIdx(pair.second)->getCompAtomType());
    mol->getAtomWithIdx(pair.first )->setCompAtomType(static_cast<int>(first_at  ^ AtomType::HbondDonor));
    mol->getAtomWithIdx(pair.second)->setCompAtomType(static_cast<int>(second_at ^ AtomType::HbondDonor));
  }

  if (mol->getRingInfo()->isInitialized()) {
    mol->getRingInfo()->reset();
  }
  mol->getRingInfo()->initialize(RDKit::FIND_RING_TYPE_SYMM_SSSR);
  mol->getRingInfo()->addAllRings(aromatic_atom_indices, aromatic_bond_indices);

}

void build_model_topology_csr(std::shared_ptr<RDKit::RWMol> &mol, RDKit::Conformer &conf, const ModelParserResult &P) {

  const auto    sequence = P.get_sequence();
  const size_t num_atoms = P.coords.size();
  const size_t num_residues = sequence.size();

  auto *info_pool = PoolFactory<InfoPool>::getFreshPoolForCurrentThread();
  auto *atom_pool = PoolFactory<AtomPool>::getFreshPoolForCurrentThread();
  auto *bond_pool = PoolFactory<BondPool>::getFreshPoolForCurrentThread();

  // total bond count to pre-allocate memory
  size_t expected_bond_count = 0;
  for (int residue_idx = 0; residue_idx < num_residues; ++residue_idx) {
    const auto &aa = sequence[residue_idx];
    const AminoAcidEdges& edges = StandardAminoAcidBondTable[aa[0]];
    expected_bond_count += edges.size + 1; // +1 for peptide bond to next residue
  }

  std::vector<RDKit::Atom*> vertices;
  vertices.reserve(num_atoms + 1);  // +1 for OXT

  std::vector<std::vector<int>> aromatic_atom_indices;
  std::vector<std::vector<int>> aromatic_bond_indices;
  std::vector<int> sulphur_atom_indices;

  // collect bonds for the CSR builder
  std::vector<RDKit::Bond*> bonds;
  bonds.reserve(expected_bond_count);

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
      const int hyb = entry.hyb[local_atom_index];
      int atom_number = StandardAminoAcidAtomicNumbers[atom_name[0]];

      auto *atom = atom_pool->createAtom();
      atom->setIdx(atom_idx);
      atom->setAtomicNum(atom_number);

      atom->setMonomerInfo(info_pool->createAtomInfo(
                atom_name,
                atom_idx + 1,
                entry.name,
                residue_idx + 1)
      );

      atom->setNumCompImplicitHs(ih);
      atom->setCompAtomType(at);
      atom->setHybridization(static_cast<RDKit::Atom::HybridizationType>(hyb));
      atom->setIdx(atom_idx);

      if (atom_number == 16 && aa[0] == 'C') { // CYS sulfur
        sulphur_atom_indices.push_back(atom_idx);
      }

      vertices.push_back(atom);
      atom_idx++;
    }

    const AminoAcidEdges& edges = StandardAminoAcidBondTable[aa[0]];

    for (size_t j = 0; j < edges.size; ++j) {
      const auto& edge = edges.edges[j];
      int atom1_idx = residue_start_idx + edge.i;
      int atom2_idx = residue_start_idx + edge.j;

      auto bond = bond_pool->createBond(atom1_idx, atom2_idx, edge.order);
      bonds.push_back(bond);
    }

    // Add peptide bond (current C to next residue N)
    if (residue_idx < static_cast<int>(num_residues) - 1) {
      int current_c_idx = residue_start_idx + 2; // C is always at position 2
      int next_n_idx    = atom_idx;              // start index of next residue (N at position 0)
      auto bond = bond_pool->createBond(current_c_idx, next_n_idx, RDKit::Bond::SINGLE);
      bonds.push_back(bond);
    }

    auto process_ring = [&aromatic_atom_indices, &vertices, residue_start_idx]
                        (const auto &atom_indices) {
      std::vector<int> ring_atom_indices;
      ring_atom_indices.reserve(atom_indices.size());
      for (int idx : atom_indices) {
        const int global_idx = residue_start_idx + idx;
        ring_atom_indices.push_back(global_idx);
        vertices[static_cast<size_t>(global_idx)]->setIsAromatic(true);
      }
      aromatic_atom_indices.push_back(std::move(ring_atom_indices));
    };

    switch (aa[0]) {
        case 'F': process_ring(phe_arom_indices); break;
        case 'Y': process_ring(tyr_arom_indices); break;
        case 'H': process_ring(his_arom_indices); break;
        case 'W':
            process_ring(trp_arom_indices5);
            process_ring(trp_arom_indices6);
            break;
    }
  }

  // Add terminal OXT atom
  auto last_entry = StandardAminoAcidDataTable[sequence.back()[0]];
  RDKit::Atom* oxt_atom = atom_pool->createAtom();
  oxt_atom->setAtomicNum(8);
  oxt_atom->setIdx(atom_idx);
  oxt_atom->setMonomerInfo(info_pool->createAtomInfo(
      "OXT",
      atom_idx + 1,
      last_entry.name,
      num_residues));
  oxt_atom->setNumCompImplicitHs(0);
  oxt_atom->setCompAtomType(9217);
  oxt_atom->setHybridization(RDKit::Atom::SP2);
  oxt_atom->setIdx(atom_idx);

  vertices.push_back(oxt_atom);

  // the first atom (N of first residue) gets 3 implicit Hs and tetrahedral geometry
  vertices.front()->setHybridization(RDKit::Atom::SP3);
  vertices.front()->setNumCompImplicitHs(3);

  // Add C-OXT bond (carbonyl carbon of last residue to OXT)
  int last_residue_c_idx = atom_idx - last_entry.size + 2;  // C is at position 2 in all AA
  {
    auto bond = bond_pool->createBond(last_residue_c_idx, atom_idx, RDKit::Bond::SINGLE);
    bonds.push_back(bond);
  }

  // for the first atom, set implicit Hs to 3
  vertices.front()->setNumCompImplicitHs(3);

  // positions
  conf.reserve(num_atoms);
  conf.setAllAtomPositions(std::move(P.coords));

  // disulfide bonds to edge list
  auto disulfide_pairs = find_disulfide_bonds(sulphur_atom_indices, conf.getPositions());
  for (const auto &pair : disulfide_pairs) {
    auto bond = bond_pool->createBond(pair.first, pair.second, RDKit::Bond::SINGLE);
    bonds.push_back(bond);
  }

  // Assign a single canonical index per logical bond before CSR build
  for (size_t i = 0; i < bonds.size(); ++i) {
    if (bonds[i]) bonds[i]->setIdx(static_cast<unsigned int>(i));
  }

  // build molecule
  mol = std::make_shared<RDKit::RWMol>(vertices, bonds, GraphType::CSRMolGraph);
  mol->addConformer(&conf, true);

  // Build ring bond indices against the final CSR graph order
  aromatic_bond_indices.clear();
  aromatic_bond_indices.reserve(aromatic_atom_indices.size());
  for (const auto &ring_atoms : aromatic_atom_indices) {
    std::vector<int> ring_bonds;
    const size_t n = ring_atoms.size();
    if (n < 3) { aromatic_bond_indices.push_back(ring_bonds); continue; }
    for (size_t i = 0; i < n; ++i) {
      int u = ring_atoms[i];
      int v = ring_atoms[(i + 1) % n];
      const RDKit::Bond *b = mol->getBondBetweenAtoms(u, v);
      if (b) {
        ring_bonds.push_back(static_cast<int>(b->getIdx()));
      } else {
        Logger::get_logger()->warn("build_model_topology_csr: No bond found between atoms {} and {} in ring", u, v);
        ring_bonds.push_back(-1);
      }
    }
    aromatic_bond_indices.push_back(std::move(ring_bonds));
  }

  for (const auto &pair : disulfide_pairs) {
    // Correct the implicit Hs from S atoms
    mol->getAtomWithIdx(pair.first) ->setNumCompImplicitHs(0);
    mol->getAtomWithIdx(pair.second)->setNumCompImplicitHs(0);

    // Correct the atom type for S-S bound atoms
    auto first_at   = static_cast<AtomType>(mol->getAtomWithIdx(pair.first) ->getCompAtomType());
    auto second_at  = static_cast<AtomType>(mol->getAtomWithIdx(pair.second)->getCompAtomType());
    mol->getAtomWithIdx(pair.first) ->setCompAtomType(static_cast<int>(first_at  ^ AtomType::HbondDonor));
    mol->getAtomWithIdx(pair.second)->setCompAtomType(static_cast<int>(second_at ^ AtomType::HbondDonor));
  }

  if (mol->getRingInfo()->isInitialized()) {
    mol->getRingInfo()->reset();
  }
  mol->getRingInfo()->initialize(RDKit::FIND_RING_TYPE_SYMM_SSSR);
  mol->getRingInfo()->addAllRings(aromatic_atom_indices, aromatic_bond_indices);
}

bool mock_build_model_topology(const ModelParserResult &P) {
  static const double MIN_COORD = -100000.0;
  static const double MAX_COORD =  100000.0;

  const auto &sequence = P.get_sequence();
  const auto &coords   = P.coords;
  if (sequence.empty() && coords.empty()) return true;

  size_t total_expected_atoms = 0;
  size_t current_coord_index  = 0;

  for (const auto &res : sequence) {
    if (res.size() != 1) { return false; }

    const char residue_letter = res[0];
    if (!StandardAminoAcidDataTable.is_valid(residue_letter)) return false;

    const auto &entry = StandardAminoAcidDataTable[residue_letter];
    size_t residue_atom_count = entry.size;

    for (size_t i = 0; i < residue_atom_count; ++i) {
      if (current_coord_index >= coords.size()) return false;

      // invalid or NaN
      const RDGeom::Point3D &pt = coords[current_coord_index];
      if (!std::isfinite(pt.x) || !std::isfinite(pt.y) || !std::isfinite(pt.z)) return false;

      // bounding box
      if (pt.x < MIN_COORD || pt.x > MAX_COORD ||
          pt.y < MIN_COORD || pt.y > MAX_COORD ||
          pt.z < MIN_COORD || pt.z > MAX_COORD) {
        return false;
      }

      ++current_coord_index;
      ++total_expected_atoms;
    }
  }

  ++total_expected_atoms; // terminal OXT (+1 atom):
  if (current_coord_index >= coords.size()) { return false; }

  const RDGeom::Point3D &oxt_pt = coords[current_coord_index];
  if (!std::isfinite(oxt_pt.x) || !std::isfinite(oxt_pt.y) || !std::isfinite(oxt_pt.z)) {
    return false;
  }
  if (oxt_pt.x < MIN_COORD || oxt_pt.x > MAX_COORD ||
      oxt_pt.y < MIN_COORD || oxt_pt.y > MAX_COORD ||
      oxt_pt.z < MIN_COORD || oxt_pt.z > MAX_COORD) {
    return false;
  }

  ++current_coord_index;
  if (total_expected_atoms != coords.size()) return false;

  return true;
}

} // namespace lahuta
