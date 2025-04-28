#include "bonds.hpp"
#include "bond_order.hpp"
#include "compute/engine.hpp"
#include "compute/result.hpp"
#include "ob/clean_mol.hpp"
#include "topology/compute.hpp"
#include "topology/data.hpp"
#include "topology/kernels.hpp"

// clang-format off
namespace lahuta::topology {

template <typename DataT>
ComputationResult BondKernel::execute(const DataContext<DataT, Mut::ReadOnly> &context, const BondComputationParams &params) {
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

    BondAssignmentResult result = assign_bonds(mol, *neighbors);

    mol.updatePropertyCache(false);
    RDKit::MolOps::setHybridization(mol);
    fix_bonds(mol);

    return ComputationResult(result);
  } catch (const std::exception &e) {
    return ComputationResult(ComputationError(std::string("Error computing bonds: ") + e.what()));
  }
}

void BondKernel::fix_bonds(RDKit::RWMol &mol) {
  for (auto atom : mol.atoms()) {
    atom->calcExplicitValence(false);
    // correct four-valent neutral N -> N+
    if (atom->getAtomicNum() == 7 && atom->getFormalCharge() == 0 && atom->getExplicitValence() == 4) {
      atom->setFormalCharge(1);
    }
  }
}

template <typename DataT>
ComputationResult NonStandardBondKernel::execute(const DataContext<DataT, Mut::ReadOnly> &context, const NonStandardBondComputationParams &params) {
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

    // NOTE: unclear if these should be free functions, but merge_bonds and fix_bonds not
    clean_bonds(result.mol, result.mol.getConformer());
    perceive_bond_orders_obabel(result.mol);
    BondKernel::fix_bonds(result.mol);

    bool include_dative_bonds = true;
    RDKit::MolOps::symmetrizeSSSR(result.mol, include_dative_bonds);
    RDKit::MolOps::setAromaticity(result.mol);

    // Merge into the original molecule
    merge_bonds(mol, result.mol, result.atom_indices);

    return ComputationResult(true);
  } catch (const std::exception &e) {
    return ComputationResult(ComputationError(std::string("Error computing non-standard bonds: ") + e.what()));
  }
}

void NonStandardBondKernel::merge_bonds(RDKit::RWMol &target, RDKit::RWMol &source, const std::vector<int> &index_map) {
  for (const auto &bond : source.bonds()) {
    int b_idx = index_map[bond->getBeginAtomIdx()];
    int e_idx = index_map[bond->getEndAtomIdx()];
    if (target.getBondBetweenAtoms(b_idx, e_idx) == nullptr) {
      auto a = target.getAtomWithIdx(b_idx);
      auto b = target.getAtomWithIdx(e_idx);
      int is_a_h = a->getAtomicNum() == 1;
      int is_b_h = b->getAtomicNum() == 1;

      if (is_a_h ^ is_b_h) {
        auto non_h_atom = a->getAtomicNum() == 1 ? b : a;
        non_h_atom->setNumExplicitHs(non_h_atom->getNumExplicitHs() + 1);
      }
      target.addBond(b_idx, e_idx, bond->getBondType());
    }
  }
}

template ComputationResult BondKernel           ::execute<TopologyData>(const DataContext<TopologyData, Mut::ReadOnly> &, const BondComputationParams &);
template ComputationResult NonStandardBondKernel::execute<TopologyData>(const DataContext<TopologyData, Mut::ReadOnly> &, const NonStandardBondComputationParams &);

} // namespace lahuta::topology
