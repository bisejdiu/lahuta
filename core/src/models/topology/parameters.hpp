#ifndef LAHUTA_MODELS_TOPOLOGY_PARAMETERS_HPP
#define LAHUTA_MODELS_TOPOLOGY_PARAMETERS_HPP

#include <GraphMol/GraphDefs.hpp>

#include "compute/parameters.hpp"

// clang-format off
namespace lahuta::models::topology {
using namespace ::lahuta::topology::compute;

namespace param_ids {
  constexpr ParameterInterface::TypeId MODEL_ATOMS        = 20;
  constexpr ParameterInterface::TypeId MODEL_BONDS        = 21;
  constexpr ParameterInterface::TypeId MODEL_POSITIONS    = 22;
  constexpr ParameterInterface::TypeId MODEL_AROMATICS    = 23;
  constexpr ParameterInterface::TypeId MODEL_DISULFIDES   = 24;
  constexpr ParameterInterface::TypeId MODEL_BUILD        = 25;
}

struct ModelAtomsParams : public ParameterBase<ModelAtomsParams> {
  static constexpr ParameterInterface::TypeId TYPE_ID = param_ids::MODEL_ATOMS;
  RDKit::GraphType graph_type = RDKit::GraphType::MolGraph;
};

struct ModelBondsParams : public ParameterBase<ModelBondsParams> {
  static constexpr ParameterInterface::TypeId TYPE_ID = param_ids::MODEL_BONDS;
};

struct ModelPositionsParams : public ParameterBase<ModelPositionsParams> {
  static constexpr ParameterInterface::TypeId TYPE_ID = param_ids::MODEL_POSITIONS;
};

struct ModelAromaticsParams : public ParameterBase<ModelAromaticsParams> {
  static constexpr ParameterInterface::TypeId TYPE_ID = param_ids::MODEL_AROMATICS;
};

struct ModelDisulfidesParams : public ParameterBase<ModelDisulfidesParams> {
  static constexpr ParameterInterface::TypeId TYPE_ID = param_ids::MODEL_DISULFIDES;
};

struct ModelBuildParams : public ParameterBase<ModelBuildParams> {
  static constexpr ParameterInterface::TypeId TYPE_ID = param_ids::MODEL_BUILD;
};

} // namespace lahuta::models::topology

#endif // LAHUTA_MODELS_TOPOLOGY_PARAMETERS_HPP
