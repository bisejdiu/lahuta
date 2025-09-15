#include "topology/engine.hpp"
#include "topology.hpp"

// clang-format off
namespace lahuta::topology {

void TopologyEngine::initialize(const TopologyBuildingOptions &opts) {
  if (auto* params = get_parameters<NeighborSearchParams>(NeighborSearchComputation<>::label)) {
    params->cutoff = opts.cutoff;
  }
  if (auto* params = get_parameters<AtomTypingParams>(AtomTypingComputation<>::label)) {
    params->use_molstar = (opts.atom_typing_method == ContactComputerType::Molstar);
  }

  enable(NonStandardBondComputation<>::label, opts.compute_nonstandard_bonds);

  engine_->set_auto_heal(opts.auto_heal);

  Logger::get_logger()->debug("TopologyEngine init: cutoff={}, auto_heal={}, nonstandard_bonds={}", opts.cutoff, opts.auto_heal, opts.compute_nonstandard_bonds);

  // TODO: we need to rethink how parameters are set and handled
  // the syntax can also be made less verbose

  // Mode gating
  if (opts.mode == TopologyBuildMode::Model) {
    enable(NeighborSearchComputation<> ::label, false);
    enable(BondComputation<>           ::label, false);
    enable(NonStandardBondComputation<>::label, false);
    enable(RingComputation<>           ::label, false);
    enable(AtomTypingComputation<>     ::label, false);

    // Virtual computation for the model-path topology
    enable(ModelTopologyComputation<>  ::label, true);
    enable(ResidueComputation<>        ::label, true);
    enable(SeedFromModelComputation<>  ::label, true);

    // Configure seeding method based on atom typing backend
    if (auto* seed_params = get_parameters<SeedFromModelParams>(SeedFromModelComputation<>::label)) {
      seed_params->use_molstar = (opts.atom_typing_method == lahuta::ContactComputerType::Molstar);
      Logger::get_logger()->debug("SeedFromModel: use_molstar={}", seed_params->use_molstar ? "true" : "false");
    }
  } else {
    // Ensure seed is off in generic mode
    enable(SeedFromModelComputation<>  ::label, false);
    enable(ModelTopologyComputation<>  ::label, false);
  }
}

} // namespace lahuta::topology
