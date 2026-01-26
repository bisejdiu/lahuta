#include <rdkit/GraphMol/RWMol.h>

#include "compute/context.hpp"
#include "compute/result.hpp"
#include "models/topology/data.hpp"
#include "models/topology/kernels.hpp"

// clang-format off
namespace lahuta::models::topology {
namespace C = lahuta::compute;

template <typename DataT>
C::ComputationResult
ModelPositionsKernel::execute(C::DataContext<DataT, C::Mut::ReadWrite> &context, const ModelPositionsParams &params) {
  auto &data = context.data();

  try {
    if (!data.conf) data.conf = std::make_unique<RDKit::Conformer>();
    data.conf->setAllAtomPositions(data.input_data->consume_coords());

    return C::ComputationResult(true);
  } catch (const std::exception &e) {
    return C::ComputationResult(C::ComputationError(std::string("Error setting model positions: ") + e.what()));
  }
}

template C::ComputationResult ModelPositionsKernel::execute<ModelData>(C::DataContext<ModelData, C::Mut::ReadWrite> &, const ModelPositionsParams &);

} // namespace lahuta::models::topology
