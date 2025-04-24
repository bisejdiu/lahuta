#include "topology.hpp"
#include "logging.hpp"

// clang-format off
namespace lahuta {

void Topology::build(TopologyBuildingOptions tops) {
  if (!mol_) {
    Logger::get_logger()->critical("Cannot build topology without a molecule.");
    throw std::runtime_error("Make sure to provide a molecule before building the topology.");
  }

  try {
    // Initialize the engine with the provided options
    engine_->initialize(tops);
    
    // Execute all computations
    bool success = engine_->execute();
    if (!success) {
    /*if (!engine_->get_engine()->run_enabled()) {*/
      Logger::get_logger()->error("Failed to execute topology computations");
      throw std::runtime_error("Failed to execute topology computations. See logs for details.");
    }
  } catch (const std::exception &e) {
    Logger::get_logger()->critical(
        "Error creating topology! Exception caught: {}. Will not terminate, "
        "but some topology-based features will be available.",
        e.what());
    throw; // Re-throw the exception to caller
  }
}

void Topology::assign_molstar_typing() {
  auto label = topology::AtomTypingComputation<>::label;
  auto* params = engine_->get_parameters<topology::AtomTypingParams>(label);
  
  if (params) {
    params->use_molstar = true;
    engine_->enable(label, true);
    engine_->execute_computation(label);
  } else {
    Logger::get_logger()->error("Could not get parameters for atom typing computation");
  }
}

void Topology::assign_arpeggio_atom_types() {
  auto label = topology::AtomTypingComputation<>::label;
  auto* params = engine_->get_parameters<topology::AtomTypingParams>(label);
  
  if (params) {
    params->use_molstar = false;
    engine_->enable(label, true);
    engine_->execute_computation(label);
  } else {
    Logger::get_logger()->error("Could not get parameters for atom typing computation");
  }
}

void Topology::enable_computation(TopologyComputation comp, bool enabled) {
  if (engine_) {
    // Handle single flag values
    if (is_base_flag(comp)) {
      engine_->enable(get_label(comp), enabled);
    } 
    // Handle combination flag values
    else {
      for (auto flag : BASE_COMPUTATION_FLAGS) {
        if (has_flag(comp, flag)) {
          engine_->enable(get_label(flag), enabled);
        }
      }
    }
  }
}

void Topology::enable_only(TopologyComputation comps) {
  if (engine_) {
    // First disable all computations
    for (auto flag : BASE_COMPUTATION_FLAGS) {
      engine_->enable(get_label(flag), false);
    }
    
    // Then enable only the requested ones
    for (auto flag : BASE_COMPUTATION_FLAGS) {
      if (has_flag(comps, flag)) {
        engine_->enable(get_label(flag), true);
      }
    }
  }
}

bool Topology::is_computation_enabled(TopologyComputation comp) const {
  if (engine_) {
    // Handle single flag values
    if (is_base_flag(comp)) {
      return engine_->is_computation_available(get_label(comp));
    }
    // Handle combination flag values
    else {
      bool result = true;
      for (auto flag : BASE_COMPUTATION_FLAGS) {
        if (has_flag(comp, flag) && !engine_->is_computation_available(get_label(flag))) {
          result = false;
        }
      }
      return result;
    }
  }
  return false;
}

bool Topology::execute_computation(TopologyComputation comp) {
  if (engine_) {
    // Handle single flag values
    if (is_base_flag(comp)) {
      return engine_->execute_computation(get_label(comp));
    }
    // Handle combination flag values
    else {
      bool success = true;
      for (auto flag : BASE_COMPUTATION_FLAGS) {
        if (has_flag(comp, flag) && !engine_->execute_computation(get_label(flag))) {
          success = false;
        }
      }
      return success;
    }
  }
  return false;
}

void Topology::set_cutoff(double cutoff) {
  if (engine_) {
    auto* params = engine_->get_parameters<topology::NeighborSearchParams>(
        topology::NeighborSearchComputation<>::label);
    if (params) {
      params->cutoff = cutoff;
    }
  }
}

void Topology::set_atom_typing_method(ContactComputerType method) {
  if (engine_) {
    auto* params = engine_->get_parameters<topology::AtomTypingParams>(
        topology::AtomTypingComputation<>::label);
    if (params) {
      params->use_molstar = (method == ContactComputerType::Molstar);
    }
  }
}

const topology::compute::ComputationLabel& Topology::get_label(TopologyComputation comp) {
  // Handle single flag values only
  switch (comp) {
    case TopologyComputation::Neighbors:
      return topology::NeighborSearchComputation<>::label;
    case TopologyComputation::Bonds:
      return topology::BondComputation<>::label;
    case TopologyComputation::Residues:
      return topology::ResidueComputation<>::label;
    case TopologyComputation::Rings:
      return topology::RingComputation<>::label;
    case TopologyComputation::AtomTyping:
      return topology::AtomTypingComputation<>::label;
    default:
      throw std::runtime_error("Invalid computation type or combination flag passed to get_label");
  }
}

size_t Topology::total_size() const {
  size_t total = sizeof(*this);

  // Get size of topology data elements
  const auto& data = engine_->get_data();
  
  /*total += sizeof(topology::AtomEntity) * data.atom_types.get_data().size();*/
  /*total += sizeof(topology::RingEntity) * data.rings.get_data().size();*/
  /*total += sizeof(topology::GroupEntity) * data.features.get_data().size();*/

  if (mol_) {
    total += sizeof(*mol_);
  }

  total += data.residues->total_size();

  return total;
}

} // namespace lahuta
