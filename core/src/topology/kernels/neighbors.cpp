#include "distances/neighbors.hpp"
#include "logging.hpp"
#include "nsgrid.hpp"
#include "topology/context.hpp"
#include "topology/kernels.hpp"

// clang-format off
namespace lahuta::topology {

template <typename DataT>
ComputationResult NeighborSearchKernel::execute(const DataContext<DataT, Mut::ReadOnly> &context, const NeighborSearchParams &params) {
  const auto &data = context.data();
  dist::NeighborSearchOptions opts;
  opts.cutoff = params.cutoff;
  auto neighbors = std::make_shared<NSResults>(dist::neighbors_within_radius_self(data.mol->getConformer().getPositions(), opts));

  Logger::get_logger()->debug("neighbors: atoms={}, cutoff={}, pairs={}", data.mol->getNumAtoms(), params.cutoff, neighbors->size());
  return ComputationResult(neighbors);
}

template ComputationResult NeighborSearchKernel::execute<TopologyContext>(const DataContext<TopologyContext, Mut::ReadOnly> &, const NeighborSearchParams &);

} // namespace lahuta::topology
