#include "topology/data.hpp"
#include "topology/kernels.hpp"
#include "aromatics.hpp"

// clang-format off
namespace lahuta::topology {

template <typename DataT>
ComputationResult
RingKernel::execute(DataContext<DataT, Mut::ReadWrite> &context, const RingComputationParams &params) {
  auto &data = context.data();
  try {
    const auto& residues = data.residues->get_residues();
    initialize_and_populate_ringinfo(*data.mol, *data.residues);
    return ComputationResult(true);
  } catch (const std::exception &e) {
    return ComputationResult(ComputationError(std::string("Error computing rings: ") + e.what()));
  }
}

template ComputationResult RingKernel::execute<TopologyData>(DataContext<TopologyData, Mut::ReadWrite> &, const RingComputationParams &);

} // namespace lahuta::topology
