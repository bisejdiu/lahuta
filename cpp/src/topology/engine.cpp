#include "topology/engine.hpp"
#include "topology.hpp"

// clang-format off
namespace lahuta::topology {

void TopologyEngine::initialize(const TopologyBuildingOptions &opts) {
  // Set parameter values based on options
  if (auto* params = get_parameters<NeighborSearchParams>(NeighborSearchComputation<>::label)) {
    params->cutoff = opts.cutoff;
  }
  
  if (auto* params = get_parameters<AtomTypingParams>(AtomTypingComputation<>::label)) {
    params->use_molstar = (opts.atom_typing_method == ContactComputerType::Molstar);
  }
  
  // Enable or disable computations based on the options
  // Currently, all computations are enabled by default
}

} // namespace lahuta::topology
