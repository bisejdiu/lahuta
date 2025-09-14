#include "aromatics.hpp"
#include "logging.hpp"
#include "topology/context.hpp"
#include "topology/kernels.hpp"

// clang-format off
namespace lahuta::topology {

template <typename DataT>
ComputationResult
RingKernel::execute(DataContext<DataT, Mut::ReadWrite> &context, const RingComputationParams &params) {
  auto &data = context.data();
  try {
    const auto& residues = data.residues->get_residues();
    initialize_and_populate_ringinfo(*data.mol, *data.residues);

    data.rings = AtomTypingKernel::populate_ring_entities(*data.mol);
    auto ring_count = data.rings.size();

    Logger::get_logger()->debug("rings: ringinfo_initialized, count={}", ring_count);
    return ComputationResult(true);
  } catch (const std::exception &e) {
    return ComputationResult(ComputationError(std::string("Error computing rings: ") + e.what()));
  }
}

template ComputationResult RingKernel::execute<TopologyContext>(DataContext<TopologyContext, Mut::ReadWrite> &, const RingComputationParams &);

} // namespace lahuta::topology
