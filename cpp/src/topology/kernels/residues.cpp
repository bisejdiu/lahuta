#include "compute/context.hpp"
#include "compute/result.hpp"
#include "topology/data.hpp"
#include "topology/kernels.hpp"
#include <rdkit/GraphMol/RWMol.h>

// clang-format off
namespace lahuta::topology {

template <typename DataT>
ComputationResult
ResidueKernel::execute(DataContext<DataT, Mut::ReadWrite> &context, const ResidueComputationParams &params) {
  auto &data = context.data();
  try {
    data.residues->build();
    return ComputationResult(true);
  } catch (const std::exception &e) {
    return ComputationResult(ComputationError(std::string("Error computing residues: ") + e.what()));
  }
}

template ComputationResult ResidueKernel::execute<TopologyData>(DataContext<TopologyData, Mut::ReadWrite> &, const ResidueComputationParams &);

} // namespace lahuta::topology
