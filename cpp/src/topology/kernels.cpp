#include "topology/kernels.hpp"
#include "aromatics.hpp"
#include "bond_order.hpp"
#include "bonds.hpp"
#include "compute/context.hpp"
#include "compute/engine.hpp"
#include "compute/result.hpp"
#include "contacts/atoms.hpp"
#include "contacts/groups.hpp"
#include "logging.hpp"
#include "ob/clean_mol.hpp"
#include "topology/compute.hpp"
#include "topology/data.hpp"
#include <rdkit/GraphMol/RWMol.h>

// clang-format off
namespace lahuta::topology {

template <typename DataT>
ComputationResult NeighborSearchKernel::execute(const DataContext<DataT, Mut::ReadOnly> &context, const NeighborSearchParams &params) {
  const auto &data = context.data();
  auto grid = FastNS(data.mol->getConformer().getPositions());
  auto ok   = grid.build(params.cutoff);
  if (!ok) return ComputationResult(ComputationError("Failed to build the grid for neighbor search."));

  auto neighbors = std::make_shared<NSResults>(grid.self_search());
  return ComputationResult(neighbors);
}

template <typename DataT>
ComputationResult BondKernel::execute(DataContext<DataT, Mut::ReadWrite> &context, const BondComputationParams &params) {
  auto &data = context.data();
  auto &mol  = *data.mol;

  try {
    // Access dependency result using the Engine
    auto engine = context.get_engine();
    if (!engine) {
      return ComputationResult(ComputationError("Engine not available for dependency resolution"));
    }

    auto maybe_neigh = engine->template result<std::shared_ptr<NSResults>>(NeighborSearchComputation<DataT>::label);
    if (!maybe_neigh)  return ComputationResult(ComputationError("Neighbor search results not available"));
    auto neighbors = *maybe_neigh;

    // Assign bonds based on neighbor results
    BondAssignmentResult result = assign_bonds(mol, *neighbors);

    // Update properties and set hybridization
    mol.updatePropertyCache(false);
    RDKit::MolOps::setHybridization(mol);
    cleanup_predef(mol);

    if (result.has_unlisted_resnames) {
      clean_bonds(result.mol, result.mol.getConformer());
      perceive_bond_orders_obabel(result.mol);
      cleanup(result.mol);

      merge_bonds(mol, result.mol, result.atom_indices);
    }

    return ComputationResult(true);
  } catch (const std::exception &e) {
    return ComputationResult(ComputationError(std::string("Error computing bonds: ") + e.what()));
  }
}

void BondKernel::cleanup_predef(RDKit::RWMol &mol) {
  for (auto atom : mol.atoms()) {
    atom->calcExplicitValence(false);
    // correct four-valent neutral N -> N+
    if (atom->getAtomicNum() == 7 && atom->getFormalCharge() == 0 && atom->getExplicitValence() == 4) {
      atom->setFormalCharge(1);
    }
  }
}

void BondKernel::cleanup(RDKit::RWMol &mol) {
  cleanup_predef(mol);

  bool include_dative_bonds = true;
  RDKit::MolOps::symmetrizeSSSR(mol, include_dative_bonds);
  RDKit::MolOps::setAromaticity(mol);
}

void BondKernel::merge_bonds(RDKit::RWMol &target, RDKit::RWMol &source, const std::vector<int> &index_map) {
  for (const auto &bond : source.bonds()) {
    int bIdx = index_map[bond->getBeginAtomIdx()];
    int eIdx = index_map[bond->getEndAtomIdx()];
    if (target.getBondBetweenAtoms(bIdx, eIdx) == nullptr) {
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

template <typename DataT>
ComputationResult
ResidueKernel::execute(DataContext<DataT, Mut::ReadWrite> &context, const ResidueComputationParams &params) {
  auto &data = context.data();
  try {
    data.residues->build();
    return ComputationResult(true);
  } catch (const std::exception &e) {
    return ComputationResult(ComputationError(std::string("Error computing residues: ") + e.what()));
  }
}

template <typename DataT>
ComputationResult
RingKernel::execute(DataContext<DataT, Mut::ReadWrite> &context, const RingComputationParams &params) {
  auto &data = context.data();
  try {
    initialize_and_populate_ringinfo(*data.mol, *data.residues);
    return ComputationResult(true);
  } catch (const std::exception &e) {
    return ComputationResult(ComputationError(std::string("Error computing rings: ") + e.what()));
  }
}

template <typename DataT>
ComputationResult
AtomTypingKernel::execute(DataContext<DataT, Mut::ReadWrite> &context, const AtomTypingParams &params) {
  auto &data = context.data();

  try {
    ValenceModel valence_model;
    valence_model.apply(*data.mol);

    data.atom_types = AtomTypeAnalysis::analyze(*data.mol);
    data.features   = GroupTypeAnalysis::analyze(*data.mol, *data.residues);
    data.rings      = populate_ring_entities(*data.mol);

    return ComputationResult(true);
  } catch (const std::exception &e) {
    Logger::get_logger()->error("Exception in atom typing: {}", e.what());
    return ComputationResult(ComputationError(std::string("Error computing atom types: ") + e.what()));
  }
}

RingEntityCollection AtomTypingKernel::populate_ring_entities(RDKit::RWMol &mol) {
  RingEntityCollection rings;
  size_t id = 0;
  for (const std::vector<int> &ring : mol.getRingInfo()->atomRings()) {
    rings.add_data(mol, ring, id++);
  }
  return rings;
}

bool AtomTypingKernel::should_initialize_ringinfo(int mol_size) {
  constexpr int small_threshold  = 20'000;
  constexpr int medium_threshold = 50'000;
  constexpr int large_threshold  = 100'000;

  if (mol_size < small_threshold) return true;

  if (mol_size < medium_threshold) {
    Logger::get_logger()->warn("Filtered molecule size ({}) is large. Performance may be affected.", mol_size);
    return true;
  } else if (mol_size < large_threshold) {
    Logger::get_logger()->warn("Filtered molecule size ({}) is very large. Performance may be severely affected.", mol_size);
    return true;
  }

  Logger::get_logger()->error("Filtered molecule size ({}) is too large. Ring perception will be skipped!", mol_size);
  return false;
}

////////// instantiate template functions
// readwrite
template ComputationResult BondKernel      ::execute<TopologyData>(DataContext<TopologyData, Mut::ReadWrite> &, const BondComputationParams &);
template ComputationResult ResidueKernel   ::execute<TopologyData>(DataContext<TopologyData, Mut::ReadWrite> &, const ResidueComputationParams &);
template ComputationResult RingKernel      ::execute<TopologyData>(DataContext<TopologyData, Mut::ReadWrite> &, const RingComputationParams &);
template ComputationResult AtomTypingKernel::execute<TopologyData>(DataContext<TopologyData, Mut::ReadWrite> &, const AtomTypingParams &);

// readonly
template ComputationResult NeighborSearchKernel::execute<TopologyData>(const DataContext<TopologyData, Mut::ReadOnly> &, const NeighborSearchParams &);

} // namespace lahuta::topology
