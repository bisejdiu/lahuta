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

    BondAssignmentResult result = assign_bonds(mol, *neighbors);

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

template ComputationResult BondKernel::execute<TopologyData>(DataContext<TopologyData, Mut::ReadWrite> &, const BondComputationParams &);

} // namespace lahuta::topology
