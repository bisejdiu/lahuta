#include <rdkit/GraphMol/MolOps.h>

#include "bonds/bonds.hpp"
#include "bonds/perception/nonstandard_bonds.hpp"
#include "bonds/perception/subset_merge.hpp"
#include "compute/engine.hpp"
#include "compute/result.hpp"
#include "logging/logging.hpp"
#include "topology/compute.hpp"
#include "topology/context.hpp"
#include "topology/kernels.hpp"

// clang-format off
namespace lahuta::topology {

template <typename DataT>
ComputationResult BondKernel::execute(DataContext<DataT, Mut::ReadWrite> &context, const BondComputationParams &params) {

  auto &data = context.data();
  auto &mol  = *data.mol;

  try {
    auto engine = context.get_engine();
    if (!engine) {
      return ComputationResult(ComputationError("Engine not available for dependency resolution"));
    }

    auto maybe_neigh = engine->template result<std::shared_ptr<NSResults>>(NeighborSearchComputation<DataT>::label);
    if (!maybe_neigh)  return ComputationResult(ComputationError("Neighbor search results not available"));
    auto neighbors = *maybe_neigh;
    Logger::get_logger()->debug("bonds: neighbors_ready, pairs={}", neighbors->size());

    BondAssignmentResult result = assign_bonds(mol, *neighbors);

    //
    // We delete neibors to free memory and avoid large memory usage. This is good from a memory standpoint, but
    // this also means that if we need to recompute neighbors for any reason later with a cutoff equal or shorter
    // than the bond cutoff, we'll have to redo the neighbor search. It's a performance vs memory tradeoff.
    // Right now we prioritize memory. - Besian, August 2025
    //
    neighbors->clear();
    neighbors.reset();

    mol.updatePropertyCache(false);
    RDKit::MolOps::setHybridization(mol);
    // fix_bonds(mol);

    Logger::get_logger()->debug("bonds: assigned, nonstandard_pending={}", result.has_unlisted_resnames);
    return ComputationResult(result);
  } catch (const std::exception &e) {
    return ComputationResult(ComputationError(std::string("Error computing bonds: ") + e.what()));
  }
}

// NOTE: this is not cheap and I doubt worth it to run calcExplicitValence
void BondKernel::fix_bonds(RDKit::RWMol &mol) {
  for (auto atom : mol.atoms()) {
    atom->calcExplicitValence(false);
    // Correct four-valent neutral N -> N+
    if (atom->getAtomicNum() == 7 && atom->getFormalCharge() == 0 && atom->getExplicitValence() == 4) {
      atom->setFormalCharge(1);
    }
  }
}

template <typename DataT>
ComputationResult NonStandardBondKernel::execute(DataContext<DataT, Mut::ReadWrite> &context, const NonStandardBondComputationParams &params) {
  auto &data = context.data();
  auto &mol  = *data.mol;

  try {
    auto engine = context.get_engine();
    if (!engine) {
      return ComputationResult(ComputationError("Engine not available for dependency resolution"));
    }

    auto maybe_result = engine->template result<BondAssignmentResult>(BondComputation<DataT>::label);
    if (!maybe_result) return ComputationResult(ComputationError("Bond computation results not available"));
    auto result = *maybe_result;

    if (!result.has_unlisted_resnames) return ComputationResult(true);

    Logger::get_logger()->debug("nonstandard_bonds: processing subset size={}", result.atom_indices.size());

    // Try residue-level cached perception
    bonds::PerceptionStats stats;
    bool used_cache = false;
    try {
      used_cache = bonds::apply_residue_level_bond_orders(result, &stats);
    } catch (const std::exception &e) {
      Logger::get_logger()->warn("nonstandard_bonds: residue-level perception failed: {}", e.what());
      used_cache = false;
    }

    if (!used_cache) {
      // Whole-subset perception. This can become expensive for large systems.
      if (!bonds::apply_whole_subset_perception(result, &stats)) {
        return ComputationResult(ComputationError("Fallback bond perception failed"));
      }
    }

    // BondKernel::fix_bonds(result.mol);

    // This should not be needed anymore, but kept temporarily in case we regress somehow
    // bool include_dative_bonds = true;
    // RDKit::MolOps::symmetrizeSSSR(result.mol, include_dative_bonds);
    // if (!used_cache) {
    //   RDKit::MolOps::setAromaticity(result.mol);
    // }

    // Merge results back into the original molecule
    bonds::subset_merge::merge_bonds(mol, result.mol, result.atom_indices);

    Logger::get_logger()->debug("nonstandard_bonds: merged");
    return ComputationResult(true);
  } catch (const std::exception &e) {
    return ComputationResult(ComputationError(std::string("Error computing non-standard bonds: ") + e.what()));
  }
}

void NonStandardBondKernel::merge_bonds(RDKit::RWMol &target, RDKit::RWMol &source, const std::vector<int> &index_map) {
  bonds::subset_merge::merge_bonds(target, source, index_map);
}

template ComputationResult BondKernel           ::execute<TopologyContext>(DataContext<TopologyContext, Mut::ReadWrite> &, const BondComputationParams &);
template ComputationResult NonStandardBondKernel::execute<TopologyContext>(DataContext<TopologyContext, Mut::ReadWrite> &, const NonStandardBondComputationParams &);

} // namespace lahuta::topology
