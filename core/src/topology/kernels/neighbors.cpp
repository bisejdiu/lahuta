#include "nsgrid.hpp"
#include "topology/data.hpp"
#include "topology/kernels.hpp"

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

template ComputationResult NeighborSearchKernel::execute<TopologyData>(const DataContext<TopologyData, Mut::ReadOnly> &, const NeighborSearchParams &);

} // namespace lahuta::topology
