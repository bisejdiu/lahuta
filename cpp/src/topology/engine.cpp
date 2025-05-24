#include "topology/engine.hpp"
#include "topology.hpp"

// clang-format off
namespace lahuta::topology {

void TopologyEngine::initialize(const TopologyBuildingOptions &opts) {
  if (auto* params = get_parameters<NeighborSearchParams>(NeighborSearchComputation<>::label)) {
    params->cutoff = opts.cutoff;
  }
  // if (auto* params = get_parameters<AtomTypingParams>(AtomTypingComputation<>::label)) {
  //   params->use_molstar = (opts.atom_typing_method == ContactComputerType::Molstar);
  // }

  enable(NonStandardBondComputation<>::label, opts.compute_nonstandard_bonds);

  engine_->set_auto_heal(opts.auto_heal);

  // TODO: we need to rethink how parameters are set and handled
  // Also, I notice, the syntax is a bit verbose.
}

} // namespace lahuta::topology
