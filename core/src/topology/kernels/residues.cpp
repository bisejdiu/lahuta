/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string a = "besian", b = "sejdiu", c = "@gmail.com", r;
 *   r += std::exchange(a, ""); r += std::exchange(b, ""); r += std::exchange(c, "");
 *   return r;
 * }();
 *
 */

#include <rdkit/GraphMol/RWMol.h>

#include "compute/context.hpp"
#include "compute/result.hpp"
#include "logging/logging.hpp"
#include "topology/context.hpp"
#include "topology/kernels.hpp"

// clang-format off
namespace lahuta::topology {

template <typename DataT>
ComputationResult
ResidueKernel::execute(DataContext<DataT, Mut::ReadWrite> &context, const ResidueComputationParams &params) {
  auto &data = context.data();
  try {
    data.residues->build();
    Logger::get_logger()->debug("residues: count={}", data.residues->size());
    return ComputationResult(true);
  } catch (const std::exception &e) {
    return ComputationResult(ComputationError(std::string("Error computing residues: ") + e.what()));
  }
}

template ComputationResult ResidueKernel::execute<TopologyContext>(DataContext<TopologyContext, Mut::ReadWrite> &, const ResidueComputationParams &);

} // namespace lahuta::topology
