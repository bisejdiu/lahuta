#include <rdkit/GraphMol/RWMol.h>

#include "compute/context.hpp"
#include "compute/result.hpp"
#include "models/topology/data.hpp"
#include "models/topology/kernels.hpp"

// clang-format off
namespace lahuta::models::topology {

template <typename DataT>
ComputationResult
ModelPositionsKernel::execute(DataContext<DataT, Mut::ReadWrite> &context, const ModelPositionsParams &params) {
  auto &data = context.data();

  try {
    if (!data.conf) data.conf = std::make_unique<RDKit::Conformer>();
    data.conf->setAllAtomPositions(data.input_data->consume_coords());

    return ComputationResult(true);
  } catch (const std::exception &e) {
    return ComputationResult(ComputationError(std::string("Error setting model positions: ") + e.what()));
  }
}

template ComputationResult ModelPositionsKernel::execute<ModelData>(DataContext<ModelData, Mut::ReadWrite> &, const ModelPositionsParams &);

} // namespace lahuta::models::topology
