#include "models/topology/compute.hpp"
#include "models/topology/engine.hpp"

// clang-format off
namespace lahuta::models::topology {

void ModelTopologyEngine::initialize(const ModelTopologyBuildingOptions &opts) {

  auto* atoms_params = get_parameters<ModelAtomsParams>(ModelAtomsComputation<>::label);
  atoms_params->graph_type = opts.graph_type;

  engine_->set_auto_heal(opts.auto_heal);

  enable_computation(ModelTopologyComputation::Atoms,      true);  // required
  enable_computation(ModelTopologyComputation::Positions,  true);  // required
  enable_computation(ModelTopologyComputation::BuildMol,   true);  // required

  enable_computation(ModelTopologyComputation::Bonds,      has_flag(opts.enabled_computations, ModelTopologyComputation::Bonds));
  enable_computation(ModelTopologyComputation::Aromatics,  has_flag(opts.enabled_computations, ModelTopologyComputation::Aromatics));
  enable_computation(ModelTopologyComputation::Disulfides, has_flag(opts.enabled_computations, ModelTopologyComputation::Disulfides));
}

void ModelTopologyEngine::enable_computation(lahuta::models::ModelTopologyComputation comp, bool enabled) {
  const auto& label = get_label(comp);
  engine_->enable(label, enabled);
}

void ModelTopologyEngine::enable_only(lahuta::models::ModelTopologyComputation comps) {
  engine_->enable(ModelAtomsComputation<>::label,       false);
  engine_->enable(ModelBondsComputation<>::label,       false);
  engine_->enable(ModelPositionsComputation<>::label,   false);
  engine_->enable(ModelAromaticsComputation<>::label,   false);
  engine_->enable(ModelDisulfidesComputation<>::label,  false);
  engine_->enable(ModelBuildComputation<>::label,       false);

  if (has_flag(comps, ModelTopologyComputation::Atoms))      { engine_->enable(ModelAtomsComputation<>::label,      true); }
  if (has_flag(comps, ModelTopologyComputation::Bonds))      { engine_->enable(ModelBondsComputation<>::label,      true); }
  if (has_flag(comps, ModelTopologyComputation::Positions))  { engine_->enable(ModelPositionsComputation<>::label,  true); }
  if (has_flag(comps, ModelTopologyComputation::Aromatics))  { engine_->enable(ModelAromaticsComputation<>::label,  true); }
  if (has_flag(comps, ModelTopologyComputation::Disulfides)) { engine_->enable(ModelDisulfidesComputation<>::label, true); }
  if (has_flag(comps, ModelTopologyComputation::BuildMol))   { engine_->enable(ModelBuildComputation<>::label,      true); }
}

bool ModelTopologyEngine::is_computation_enabled(lahuta::models::ModelTopologyComputation comp) const {
  const auto& label = get_label(comp);
  return is_computation_available(label);
}

const ComputationLabel& ModelTopologyEngine::get_label(lahuta::models::ModelTopologyComputation comp) {
  switch (comp) {
    case ModelTopologyComputation::Atoms:      return ModelAtomsComputation<>::label;
    case ModelTopologyComputation::Bonds:      return ModelBondsComputation<>::label;
    case ModelTopologyComputation::Positions:  return ModelPositionsComputation<>::label;
    case ModelTopologyComputation::Aromatics:  return ModelAromaticsComputation<>::label;
    case ModelTopologyComputation::Disulfides: return ModelDisulfidesComputation<>::label;
    case ModelTopologyComputation::BuildMol:   return ModelBuildComputation<>::label;
    default:
      throw std::runtime_error("Invalid ModelTopologyComputation");
  }
}

} // namespace lahuta::models::topology
