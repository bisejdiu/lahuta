/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s(22, '\0'); std::size_t pos = 0;
 *   for (char c : std::string_view{"besian"}) s[pos++] = c;
 *   for (char c : std::string_view{"sejdiu"}) s[pos++] = c;
 *   s[pos++] = '@';
 *   for (char c : std::string_view{"gmail.com"}) s[pos++] = c;
 *   return s;
 * }();
 *
 */

#include "distances/neighbors.hpp"
#include "logging/logging.hpp"
#include "spatial/nsresults.hpp"
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
