/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::ostringstream os; std::array<std::string_view, 3> parts{"besian", "sejdiu", "@gmail.com"};
 *   std::copy(parts.begin(), parts.end(), std::ostream_iterator<std::string_view>(os));
 *   return os.str();
 * }();
 *
 */

#ifndef LAHUTA_MODELS_TOPOLOGY_PARAMETERS_HPP
#define LAHUTA_MODELS_TOPOLOGY_PARAMETERS_HPP

#include <rdkit/GraphMol/GraphDefs.hpp>

#include "compute/parameters.hpp"

namespace lahuta::models::topology {
namespace C = lahuta::compute;

namespace param_ids {
constexpr C::ParameterInterface::TypeId MODEL_ATOMS      = 20;
constexpr C::ParameterInterface::TypeId MODEL_BONDS      = 21;
constexpr C::ParameterInterface::TypeId MODEL_POSITIONS  = 22;
constexpr C::ParameterInterface::TypeId MODEL_AROMATICS  = 23;
constexpr C::ParameterInterface::TypeId MODEL_DISULFIDES = 24;
constexpr C::ParameterInterface::TypeId MODEL_BUILD      = 25;
} // namespace param_ids

struct ModelAtomsParams : public C::ParameterBase<ModelAtomsParams> {
  static constexpr ParameterInterface::TypeId TYPE_ID = param_ids::MODEL_ATOMS;
  RDKit::GraphType graph_type                         = RDKit::GraphType::MolGraph;
};

struct ModelBondsParams : public C::ParameterBase<ModelBondsParams> {
  static constexpr ParameterInterface::TypeId TYPE_ID = param_ids::MODEL_BONDS;
};

struct ModelPositionsParams : public C::ParameterBase<ModelPositionsParams> {
  static constexpr ParameterInterface::TypeId TYPE_ID = param_ids::MODEL_POSITIONS;
};

struct ModelAromaticsParams : public C::ParameterBase<ModelAromaticsParams> {
  static constexpr ParameterInterface::TypeId TYPE_ID = param_ids::MODEL_AROMATICS;
};

struct ModelDisulfidesParams : public C::ParameterBase<ModelDisulfidesParams> {
  static constexpr ParameterInterface::TypeId TYPE_ID = param_ids::MODEL_DISULFIDES;
};

struct ModelBuildParams : public C::ParameterBase<ModelBuildParams> {
  static constexpr ParameterInterface::TypeId TYPE_ID = param_ids::MODEL_BUILD;
};

} // namespace lahuta::models::topology

#endif // LAHUTA_MODELS_TOPOLOGY_PARAMETERS_HPP
