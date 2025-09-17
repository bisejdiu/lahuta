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
    Logger::get_logger()->debug("Topology build start");

    engine_->initialize(tops);

    if (tops.mode == TopologyBuildMode::Model) {
      bool conflicts = false;
      conflicts |= engine_->is_computation_available(topology::NeighborSearchComputation<>::label);
      conflicts |= engine_->is_computation_available(topology::BondComputation<>::label);
      conflicts |= engine_->is_computation_available(topology::NonStandardBondComputation<>::label);
      conflicts |= engine_->is_computation_available(topology::RingComputation<>::label);
      conflicts |= engine_->is_computation_available(topology::AtomTypingComputation<>::label);
      if (conflicts) {
        Logger::get_logger()->error("Model mode requires only residues + seed_from_model to be enabled");
        throw std::runtime_error("Invalid computation configuration for Model mode");
      }
    }

    bool success = engine_->execute();
    if (!success) {
      Logger::get_logger()->error("Failed to execute topology computations");
      throw std::runtime_error("Failed to execute topology computations. See logs for details.");
    }
    Logger::get_logger()->debug("Topology build done");
  } catch (const std::exception &e) {
    Logger::get_logger()->critical(
        "Error creating topology! Exception caught: {} Will not terminate, "
        "but some topology-based features will not be available.", e.what());
    throw;
  }
}

void Topology::assign_typing(AtomTypingMethod method) {
  auto& label  = topology::AtomTypingComputation<>::label;
  auto* params = engine_->get_parameters<topology::AtomTypingParams>(label);

  if (params) {
    params->mode = method;
    engine_->enable(label, true);

    Logger::get_logger()->debug("Executing atom typing ({})", contact_computer_name(method));
    engine_->execute_computation(label);

  } else {
    Logger::get_logger()->error("Could not get parameters for atom typing computation");
  }
}

void Topology::enable_computation(TopologyComputation comp, bool enabled) {
  if (!engine_) throw std::runtime_error("No engine available");

  if (is_base_flag(comp)) {
    Logger::get_logger()->debug("{} computation {}", get_label(comp).to_string_view(), enabled ? "enabled" : "disabled");
    engine_->enable(get_label(comp), enabled);
    return;
  }

  for (auto flag : BASE_COMPUTATION_FLAGS) {
    if (has_flag(comp, flag)) {
      Logger::get_logger()->debug("{} computation {}", get_label(flag).to_string_view(), enabled ? "enabled" : "disabled");
      engine_->enable(get_label(flag), enabled);
    }
  }
}

void Topology::enable_only(TopologyComputation comps) {
  if (!engine_) throw std::runtime_error("No engine available");

  for (auto flag : BASE_COMPUTATION_FLAGS) {
    engine_->enable(get_label(flag), false);
  }

  for (auto flag : BASE_COMPUTATION_FLAGS) {
    if (has_flag(comps, flag)) {
      engine_->enable(get_label(flag), true);
    }
  }
}

bool Topology::is_computation_enabled(TopologyComputation comp) const {
  if (!engine_) throw std::runtime_error("No engine available");
  if (is_base_flag(comp)) return engine_->is_computation_available(get_label(comp));

  bool result = true;
  for (auto flag : BASE_COMPUTATION_FLAGS) {
    if (has_flag(comp, flag) && !engine_->is_computation_available(get_label(flag))) {
      result = false;
    }
  }
  return result;
}

bool Topology::execute_computation(TopologyComputation comp) {
  if (!engine_) return false;

  if (is_base_flag(comp)) {
    Logger::get_logger()->debug("Running computation: {}", get_label(comp).to_string_view());
    return engine_->execute_computation(get_label(comp));
  }

  bool success = true;
  for (auto flag : BASE_COMPUTATION_FLAGS) {
    if (has_flag(comp, flag)) {
      Logger::get_logger()->debug("Running computation: {}", get_label(flag).to_string_view());
      success &= engine_->execute_computation(get_label(flag));
    }
  }

  return success;
}

void Topology::set_cutoff(double cutoff) {
  if (!engine_) throw std::runtime_error("No engine available");
  auto* params = engine_->get_parameters<topology::NeighborSearchParams>(topology::NeighborSearchComputation<>::label);
  if (params) params->cutoff = cutoff;
}

void Topology::set_atom_typing_method(AtomTypingMethod method) {
  if (!engine_) throw std::runtime_error("No engine available");
  auto* params = engine_->get_parameters<topology::AtomTypingParams>(topology::AtomTypingComputation<>::label);
  if (params) params->mode = method;
}

void Topology::set_compute_nonstandard_bonds(bool compute) {
  if (!engine_) throw std::runtime_error("No engine available");
  engine_->enable(topology::NonStandardBondComputation<>::label, compute);
}

const topology::ComputationLabel& Topology::get_label(TopologyComputation comp) {
  switch (comp) {
    // this is a bit hidden here and we may easily forget to update it if we add new computations
    case TopologyComputation::Neighbors:         return topology::NeighborSearchComputation<>::label;
    case TopologyComputation::Bonds:             return topology::BondComputation<>::label;
    case TopologyComputation::NonStandardBonds:  return topology::NonStandardBondComputation<>::label;
    case TopologyComputation::Residues:          return topology::ResidueComputation<>::label;
    case TopologyComputation::Rings:             return topology::RingComputation<>::label;
    case TopologyComputation::AtomTyping:        return topology::AtomTypingComputation<>::label;
    default:
      throw std::runtime_error("Invalid computation type or combination flag passed to get_label");
  }
}

// FIX: We should update our code so we do not rely on this anymore.
size_t Topology::total_size() const {
  size_t total = sizeof(*this);

  // Get size of topology data elements
  const auto& data = engine_->get_data();

  total += sizeof(AtomRec)  * data.atoms.size();
  total += sizeof(RingRec)  * data.rings.size();
  total += sizeof(GroupRec) * data.groups.size();

  total += sizeof(AtomRec)  * engine_->get_data().atoms.size();
  total += sizeof(RingRec)  * engine_->get_data().rings.size();
  total += sizeof(GroupRec) * engine_->get_data().groups.size();

  if (mol_) { total += sizeof(*mol_); }
  total += data.residues->total_size();

  return total;
}

const AtomRec& Topology::atom(uint32_t idx) const {
  return resolve<Kind::Atom>(EntityID::make(Kind::Atom, idx));
}

const RingRec& Topology::ring(uint32_t idx) const {
  return resolve<Kind::Ring>(EntityID::make(Kind::Ring, idx));
}

const GroupRec& Topology::group(uint32_t idx) const {
  return resolve<Kind::Group>(EntityID::make(Kind::Group, idx));
}

} // namespace lahuta
