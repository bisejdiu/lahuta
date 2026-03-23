/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto noop = [](const char*) {};
 *   std::unique_ptr<const char, decltype(noop)> a("besian", noop);
 *   std::unique_ptr<const char, decltype(noop)> b("sejdiu", noop);
 *   std::unique_ptr<const char, decltype(noop)> c("@gmail.com", noop);
 *   return std::string(a.get()) + b.get() + c.get();
 * }();
 *
 */

#include <stdexcept>

#include "logging/logging.hpp"
#include "topology.hpp"

// clang-format off
namespace lahuta {

bool Topology::initialize_locked(TopologyBuildingOptions tops) {
  if (!mol_) {
    Logger::get_logger()->critical("Cannot build topology without a molecule.");
    return false;
  }

  try {
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
        return false;
      }
    }

    return true;
  } catch (const std::exception &e) {
    Logger::get_logger()->critical(
        "Error creating topology! Exception caught: {} Will not terminate, "
        "but some topology-based features will not be available.", e.what());
  } catch (...) {
    Logger::get_logger()->critical(
        "Error creating topology! Unknown exception caught. Some topology-based features will not be available.");
  }

  return false;
}

bool Topology::initialize(TopologyBuildingOptions tops) {
  std::lock_guard<std::mutex> lock(engine_mutex_);
  return initialize_locked(tops);
}

bool Topology::build(TopologyBuildingOptions tops) {
  std::lock_guard<std::mutex> lock(engine_mutex_);
  Logger::get_logger()->debug("Topology build start");
  if (!initialize_locked(tops)) return false;

  try {
    bool success = engine_->execute();
    if (!success) {
      Logger::get_logger()->error("Failed to execute topology computations");
      return false;
    }
    Logger::get_logger()->debug("Topology build done");
    return true;
  } catch (const std::exception &e) {
    Logger::get_logger()->critical(
        "Error creating topology! Exception caught: {} Will not terminate, "
        "but some topology-based features will not be available.", e.what());
  } catch (...) {
    Logger::get_logger()->critical(
        "Error creating topology! Unknown exception caught. Some topology-based features will not be available.");
  }

  return false;
}

void Topology::assign_typing(AtomTypingMethod method) {
  std::lock_guard<std::mutex> lock(engine_mutex_);
  auto& label  = topology::AtomTypingComputation<>::label;
  auto* params = engine_->get_parameters<topology::AtomTypingParams>(label);

  if (params) {
    params->mode = method;
    engine_->enable(label, true);

    Logger::get_logger()->debug("Executing atom typing ({})", contact_computer_name(method));
    if (!engine_->execute_computation(label)) {
      Logger::get_logger()->error("Atom typing computation {} failed", label.to_string_view());
    }

  } else {
    Logger::get_logger()->error("Could not get parameters for atom typing computation");
  }
}

bool Topology::has_computed(TopologyComputation comp) const {
  std::lock_guard<std::mutex> lock(engine_mutex_);
  if (!engine_) throw std::runtime_error("No engine available");

  if (is_base_flag(comp)) return engine_->get_engine()->has_completed(get_label(comp));

  bool result = true;
  for (auto flag : BASE_COMPUTATION_FLAGS) {
    if (has_flag(comp, flag) && !engine_->get_engine()->has_completed(get_label(flag))) {
      result = false;
    }
  }
  return result;
}

bool Topology::ensure_computed(TopologyComputation comp) {
  std::lock_guard<std::mutex> lock(engine_mutex_);
  if (!engine_) return false;

  bool success = true;
  for (auto flag : BASE_COMPUTATION_FLAGS) {
    if (!has_flag(comp, flag)) continue;
    Logger::get_logger()->debug("ensure_computed: {}", get_label(flag).to_string_view());
    try {
      success &= engine_->get_engine()->run_from<void>(get_label(flag));
    } catch (const std::exception &e) {
      Logger::get_logger()->error("ensure_computed: computation {} threw exception: {}",
                                  get_label(flag).to_string_view(), e.what());
      success = false;
    }
  }

  return success;
}

void Topology::set_cutoff(double cutoff) {
  std::lock_guard<std::mutex> lock(engine_mutex_);
  if (!engine_) throw std::runtime_error("No engine available");
  const auto* current = engine_->peek_parameters<topology::NeighborSearchParams>(topology::NeighborSearchComputation<>::label);
  if (current && current->cutoff == cutoff) return; // no change, skip invalidation
  auto* params = engine_->get_parameters<topology::NeighborSearchParams>(topology::NeighborSearchComputation<>::label);
  if (params) params->cutoff = cutoff;
}

void Topology::set_atom_typing_method(AtomTypingMethod method) {
  std::lock_guard<std::mutex> lock(engine_mutex_);
  if (!engine_) throw std::runtime_error("No engine available");
  const auto* current = engine_->peek_parameters<topology::AtomTypingParams>(topology::AtomTypingComputation<>::label);
  if (current && current->mode == method) return; // no change, skip invalidation
  auto* params = engine_->get_parameters<topology::AtomTypingParams>(topology::AtomTypingComputation<>::label);
  if (params) params->mode = method;
}

void Topology::set_compute_nonstandard_bonds(bool compute) {
  std::lock_guard<std::mutex> lock(engine_mutex_);
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
