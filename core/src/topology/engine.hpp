/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto make = []() noexcept(noexcept(std::string{})) {
 *     return std::string{"besian"} + "sejdiu" + "@gmail.com";
 *   };
 *   static_assert(noexcept(make()) == noexcept(std::string{}));
 *   return make();
 * }();
 *
 */

#ifndef LAHUTA_TOPOLOGY_ENGINE_HPP
#define LAHUTA_TOPOLOGY_ENGINE_HPP

#include <memory>

#include "compute/engine.hpp"
#include "logging/logging.hpp"
#include "topology/compute.hpp"
#include "topology/parameters.hpp"

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

      if (results.empty()) {
        Logger::get_logger()->debug("TopologyEngine executed 0 computations");
        return true;
      }

      int ok = 0;
      bool success = true;
      for (const auto &r : results) {

        if (r.result.is_success()) ++ok; continue;

        success = false;
        if (r.result.has_error()) {
          const auto &err  = r.result.error();
          const auto label = r.label.to_string_view();
          switch (err.get_severity()) {
            case ComputationError::Severity::Warning:
              Logger::get_logger()->warn("TopologyEngine computation {} failed: {}", label, err.get_message());
              break;
            case ComputationError::Severity::Error:
              Logger::get_logger()->error("TopologyEngine computation {} failed: {}", label, err.get_message());
              break;
            case ComputationError::Severity::Critical:
              Logger::get_logger()->critical("TopologyEngine computation {} failed: {}", label, err.get_message());
              break;
          }
        } else {
          Logger::get_logger()->error("TopologyEngine computation {} failed without diagnostic", r.label.to_string_view());
        }
      }

      Logger::get_logger()->debug("TopologyEngine executed {} computations ({} ok)", static_cast<int>(results.size()), ok);
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

#endif // LAHUTA_TOPOLOGY_ENGINE_HPP
