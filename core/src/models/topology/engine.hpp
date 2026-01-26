#ifndef LAHUTA_MODELS_TOPOLOGY_ENGINE_HPP
#define LAHUTA_MODELS_TOPOLOGY_ENGINE_HPP

#include <algorithm>
#include <memory>

#include "compute/engine.hpp"
#include "logging/logging.hpp"
#include "models/topology/compute.hpp"
#include "models/topology/data.hpp"
#include "models/topology/parameters.hpp"

namespace lahuta::models {

enum class ModelTopologyComputation : uint8_t { // TODO: needs more testing
  None       = 0,
  Atoms      = 1 << 0,
  Bonds      = 1 << 1,
  Positions  = 1 << 2,
  Aromatics  = 1 << 3,
  Disulfides = 1 << 4,
  BuildMol   = 1 << 5,
  All        = Atoms | Bonds | Positions | Aromatics | Disulfides | BuildMol
};

inline ModelTopologyComputation operator|(ModelTopologyComputation a, ModelTopologyComputation b) {
  return static_cast<ModelTopologyComputation>(static_cast<uint8_t>(a) | static_cast<uint8_t>(b));
}

inline bool has_flag(ModelTopologyComputation flags, ModelTopologyComputation flag) {
  return (static_cast<uint8_t>(flags) & static_cast<uint8_t>(flag)) != 0;
}

struct ModelTopologyBuildingOptions {
  RDKit::GraphType graph_type = RDKit::GraphType::CSRMolGraph;
  bool auto_heal              = true;

  ModelTopologyComputation enabled_computations = ModelTopologyComputation::All;
};

} // namespace lahuta::models

namespace lahuta::models::topology {

// clang-format off
class ModelTopologyEngine {
public:
  explicit ModelTopologyEngine(const ModelParserResult& input) {
    data_   = std::make_unique<ModelData>(input);
    engine_ = std::make_unique<C::ComputeEngine<ModelData, C::Mut::ReadWrite>>(*data_);

    register_computations();

    engine_->enable(ModelAtomsComputation<>::label,        true);  // required
    engine_->enable(ModelPositionsComputation<>::label,    true);  // required
    engine_->enable(ModelBuildComputation<>::label,        true);  // required
    engine_->enable(ModelBondsComputation<>::label,        true);  // optional
    engine_->enable(ModelAromaticsComputation<>::label,    true);  // optional
    engine_->enable(ModelDisulfidesComputation<>::label,   true);  // optional
  }

  void initialize(const ModelTopologyBuildingOptions &opts);

  /// Execute all enabled computations
  bool execute() {
    try {
      auto results    = engine_->w_execute_all();
      auto is_success = [](const C::ResultEntry &entry) { return entry.result.is_success(); };
      bool success    = !results.empty() && std::all_of(results.begin(), results.end(), is_success);

      return success;
    } catch (const std::exception &e) {
      Logger::get_logger()->error("Error executing model topology computations: {}", e.what());
      return false;
    }
  }

  bool execute_computation(const C::ComputationLabel &label) {
    try {
      return engine_->run<void>(label);
    } catch (const std::exception &e) {
      Logger::get_logger()->error("Error executing model computation {}: {}",
                                  label.to_string_view(),
                                  e.what());
      return false;
    }
  }

  ModelData &get_data() { return *data_; }
  const ModelData &get_data() const { return *data_; }

  void enable(const C::ComputationLabel &label, bool enabled) { engine_->enable(label, enabled); }

  void enable_computation(ModelTopologyComputation comp, bool enabled);
  void enable_only(ModelTopologyComputation comps);
  bool is_computation_enabled(ModelTopologyComputation comp) const;

  bool is_computation_available(const C::ComputationLabel &label) const {
    int idx = engine_->find_label(label);
    return idx >= 0 && engine_->is_enabled(idx);
  }

  template <typename P>
  P *get_parameters(const C::ComputationLabel &label) {
    return &(engine_->template get_parameters<P>(label));
  }

  C::ComputeEngine<ModelData, C::Mut::ReadWrite> *get_engine() {
    return static_cast<C::ComputeEngine<ModelData, C::Mut::ReadWrite> *>(engine_.get());
  }

private:
  void register_computations() {
    engine_->add(std::make_unique<ModelAtomsComputation<>>(ModelAtomsParams{}));
    engine_->add(std::make_unique<ModelBondsComputation<>>(ModelBondsParams{}));
    engine_->add(std::make_unique<ModelPositionsComputation<>>(ModelPositionsParams{}));
    engine_->add(std::make_unique<ModelAromaticsComputation<>>(ModelAromaticsParams{}));
    engine_->add(std::make_unique<ModelDisulfidesComputation<>>(ModelDisulfidesParams{}));
    engine_->add(std::make_unique<ModelBuildComputation<>>(ModelBuildParams{}));
  }

  static const C::ComputationLabel &get_label(ModelTopologyComputation comp);

private:
  std::unique_ptr<ModelData> data_;
  std::unique_ptr<C::ComputeEngine<ModelData, C::Mut::ReadWrite>> engine_;
};

} // namespace lahuta::models::topology

#endif // LAHUTA_MODELS_TOPOLOGY_ENGINE_HPP
