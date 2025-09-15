#pragma once

#include "compute/engine.hpp"
#include "logging.hpp"
#include "topology/compute.hpp"
#include "topology/parameters.hpp"
#include <memory>

// clang-format off
namespace lahuta {
struct TopologyBuildingOptions;
} // namespace lahuta

namespace lahuta::topology {
using namespace compute;

/// manages the computation engine for topology
class TopologyEngine {
public:
  explicit TopologyEngine(std::shared_ptr<RDKit::RWMol> mol) {
    data_   = std::make_unique<TopologyContext>(mol);
    engine_ = std::make_unique<ComputeEngine<TopologyContext, Mut::ReadWrite>>(*data_);

    register_computation();

    engine_->set_auto_heal(true); // TODO: replace with TopologyBuildingOptions parameter

    engine_->enable(NeighborSearchComputation<> ::label, true);
    engine_->enable(BondComputation<>           ::label, true);
    engine_->enable(NonStandardBondComputation<>::label, true);
    engine_->enable(ResidueComputation<>        ::label, true);
    engine_->enable(RingComputation<>           ::label, true);
    engine_->enable(AtomTypingComputation<>     ::label, true);
    engine_->enable(SeedFromModelComputation<>  ::label, false);
  }

  void initialize(const TopologyBuildingOptions &opts);

  /// Execute all enabled computations
  bool execute() {
    try {
      auto results = engine_->w_execute_all();

      // code block is just for logging
      if (!results.empty()) {
        int ok = 0;
        for (const auto &r : results) ok += r.result.is_success() ? 1 : 0;
        Logger::get_logger()->debug("TopologyEngine executed {} computations ({} ok)", (int)results.size(), ok);
        for (const auto &r : results) {
          if (!r.result.is_success()) {
            Logger::get_logger()->debug("  - {}: failed", r.label.to_string_view());
          }
        }
      } else {
        Logger::get_logger()->debug("TopologyEngine executed 0 computations");
      }
      auto is_success = [](const ResultEntry &entry) { return entry.result.is_success(); };
      bool success = !results.empty() && std::all_of(results.begin(), results.end(), is_success);

      return success;
    } catch (const std::exception &e) {
      Logger::get_logger()->error("Error executing topology computations: {}", e.what());
      return false;
    }
  }

  bool execute_computation(const ComputationLabel &label) {
    try {
      Logger::get_logger()->debug("TopologyEngine run: {}", label.to_string_view());
      return engine_->run<void>(label);
    } catch (const std::exception &e) {
      Logger::get_logger()->error("Error executing computation {}: {}", label.to_string_view(), e.what());
      return false;
    }
  }

  TopologyContext &get_data() { return *data_; }
  const TopologyContext &get_data() const { return *data_; }

  void enable(const ComputationLabel &label, bool enabled) {
    engine_->enable(label, enabled);
  }

  bool is_computation_available(const ComputationLabel &label) const {
    int idx = engine_->find_label(label);
    return idx >= 0 && engine_->is_enabled(idx);
  }

  template <typename P>
  P *get_parameters(const ComputationLabel &label) {
    return &(engine_->template get_parameters<P>(label));
  }

  ComputeEngine<TopologyContext, Mut::ReadWrite> *get_engine() {
    return static_cast<ComputeEngine<TopologyContext, Mut::ReadWrite> *>(engine_.get());
  }

private:
  void register_computation() {
    engine_->add(std::make_unique<NeighborSearchComputation<>>(NeighborSearchParams{}));
    engine_->add(std::make_unique<BondComputation<>>(BondComputationParams{}));
    engine_->add(std::make_unique<NonStandardBondComputation<>>(NonStandardBondComputationParams{}));
    engine_->add(std::make_unique<ResidueComputation<>>(ResidueComputationParams{}));
    engine_->add(std::make_unique<RingComputation<>>(RingComputationParams{}));
    engine_->add(std::make_unique<AtomTypingComputation<>>(AtomTypingParams{}));
    engine_->add(std::make_unique<SeedFromModelComputation<>>(SeedFromModelParams{}));
    engine_->add(std::make_unique<ModelTopologyComputation<>>(ModelTopologyParams{}));
  }

private:
  std::unique_ptr<TopologyContext> data_;
  std::unique_ptr<ComputeEngine<TopologyContext, Mut::ReadWrite>> engine_;
};

} // namespace lahuta::topology
