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
