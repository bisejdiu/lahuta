#pragma once

#include "compute/context.hpp"
#include "compute/result.hpp"
#include "entities.hpp"
#include "parameters.hpp"
#include <rdkit/GraphMol/RWMol.h>

namespace lahuta::topology {
using namespace compute;

struct NeighborSearchKernel {
  template <typename DataT>
  static ComputationResult
  execute(const DataContext<DataT, Mut::ReadOnly> &context, const NeighborSearchParams &params);
};

struct BondKernel {
  template <typename DataT>
  static ComputationResult
  execute(DataContext<DataT, Mut::ReadWrite> &context, const BondComputationParams &params);

private:
  static void cleanup_predef(RDKit::RWMol &mol);
  static void cleanup(RDKit::RWMol &mol);
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

private:
  static RingEntityCollection populate_ring_entities(RDKit::RWMol &mol);
  static bool should_initialize_ringinfo(int mol_size);
};

} // namespace lahuta::topology
