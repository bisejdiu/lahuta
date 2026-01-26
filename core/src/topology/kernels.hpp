#ifndef LAHUTA_TOPOLOGY_KERNELS_HPP
#define LAHUTA_TOPOLOGY_KERNELS_HPP

#include <rdkit/GraphMol/RWMol.h>

#include "compute/context.hpp"
#include "compute/result.hpp"
#include "entities/records.hpp"
#include "parameters.hpp"

namespace lahuta::topology {
using namespace compute;

// For virtual/aggregate computations whose execution is represented purely by dependencies in the compute graph.
struct NoopKernel {
  template <typename DataT, typename ParamT>
  static ComputationResult
  execute(const DataContext<DataT, Mut::ReadOnly> &/*context*/, const ParamT &/*params*/) {
    return ComputationResult(true);
  }
};

struct NeighborSearchKernel {
  template <typename DataT>
  static ComputationResult
  execute(const DataContext<DataT, Mut::ReadOnly> &context, const NeighborSearchParams &params);
};

struct BondKernel {
  template <typename DataT>
  static ComputationResult
  execute(DataContext<DataT, Mut::ReadWrite> &context, const BondComputationParams &params);
  static void fix_bonds(RDKit::RWMol &mol);
};

struct NonStandardBondKernel {
  template <typename DataT>
  static ComputationResult
  execute(DataContext<DataT, Mut::ReadWrite> &context, const NonStandardBondComputationParams &params);

private:
  static void merge_bonds(RDKit::RWMol &target, RDKit::RWMol &source, const std::vector<int> &index_map);
};

struct ResidueKernel {
  template <typename DataT>
  static ComputationResult
  execute(DataContext<DataT, Mut::ReadWrite> &context, const ResidueComputationParams &params);
};

struct RingKernel {
  template <typename DataT>
  static ComputationResult
  execute(DataContext<DataT, Mut::ReadWrite> &context, const RingComputationParams &params);
};

struct AtomTypingKernel {
  template <typename DataT>
  static ComputationResult
  execute(DataContext<DataT, Mut::ReadWrite> &context, const AtomTypingParams &params);

  static std::vector<RingRec> populate_ring_entities(RDKit::RWMol &mol);
private:
  static bool should_initialize_ringinfo(int mol_size);
};

struct SeedFromModelKernel {
  template <typename DataT>
  static ComputationResult
  execute(DataContext<DataT, Mut::ReadWrite> &context, const SeedFromModelParams &params);
};

} // namespace lahuta::topology

#endif // LAHUTA_TOPOLOGY_KERNELS_HPP
