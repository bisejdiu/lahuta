#include "logging.hpp"
#include "models/topology/compute.hpp"
#include "models/topology_impl.hpp"

// clang-format off
namespace lahuta::models {

bool ModelTopology::build(const ModelTopologyBuildingOptions& options) {
  engine_->initialize(options);

  const bool success = engine_->execute();
  if (!success) {
    Logger::get_logger()->error("Failed to build model topology");
  }

  return success;
}

void ModelTopology::run_computations(ModelTopologyComputation mask) {
  using namespace topology;

  if (has_flag(mask, ModelTopologyComputation::Atoms))      engine_->execute_computation(ModelAtomsComputation<>::label);
  if (has_flag(mask, ModelTopologyComputation::Bonds))      engine_->execute_computation(ModelBondsComputation<>::label);
  if (has_flag(mask, ModelTopologyComputation::Positions))  engine_->execute_computation(ModelPositionsComputation<>::label);
  if (has_flag(mask, ModelTopologyComputation::Aromatics))  engine_->execute_computation(ModelAromaticsComputation<>::label);
  if (has_flag(mask, ModelTopologyComputation::Disulfides)) engine_->execute_computation(ModelDisulfidesComputation<>::label);
  if (has_flag(mask, ModelTopologyComputation::BuildMol))   engine_->execute_computation(ModelBuildComputation<>::label);
}

void ModelTopology::enable_computation(ModelTopologyComputation comp, bool enabled) {
  engine_->enable_computation(comp, enabled);
}

void ModelTopology::enable_only(ModelTopologyComputation comps) {
  engine_->enable_only(comps);
}

bool ModelTopology::is_computation_enabled(ModelTopologyComputation comp) const {
  return engine_->is_computation_enabled(comp);
}

bool ModelTopology::execute_computation(ModelTopologyComputation comp) {
  using namespace topology;

  switch (comp) {
    case ModelTopologyComputation::Atoms:      return engine_->execute_computation(ModelAtomsComputation<>::label);
    case ModelTopologyComputation::Bonds:      return engine_->execute_computation(ModelBondsComputation<>::label);
    case ModelTopologyComputation::Positions:  return engine_->execute_computation(ModelPositionsComputation<>::label);
    case ModelTopologyComputation::Aromatics:  return engine_->execute_computation(ModelAromaticsComputation<>::label);
    case ModelTopologyComputation::Disulfides: return engine_->execute_computation(ModelDisulfidesComputation<>::label);
    case ModelTopologyComputation::BuildMol:   return engine_->execute_computation(ModelBuildComputation<>::label);
    default:
      return false;
  }
}

void ModelTopology::set_graph_type(RDKit::GraphType type) {
  auto* params = engine_->get_parameters<topology::ModelAtomsParams>(topology::ModelAtomsComputation<>::label);
  params->graph_type = type;
}

} // namespace lahuta::models
