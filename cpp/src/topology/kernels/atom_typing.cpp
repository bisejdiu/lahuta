#include "compute/context.hpp"
#include "compute/result.hpp"
#include "contacts/atoms.hpp"
#include "contacts/groups.hpp"
#include "logging.hpp"
#include "topology/data.hpp"
#include "topology/kernels.hpp"
#include <rdkit/GraphMol/RWMol.h>
#include <valence_model.hpp>

// clang-format off
namespace lahuta::topology {

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

template ComputationResult AtomTypingKernel::execute<TopologyData>(DataContext<TopologyData, Mut::ReadWrite> &, const AtomTypingParams &);

} // namespace lahuta::topology
