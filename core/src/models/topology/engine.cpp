/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: namespace detail_c46 {
 *   constexpr std::array<const char*, 3> parts{"besian", "sejdiu", "@gmail.com"};
 *   template<std::size_t... Is>
 *   std::string expand(std::index_sequence<Is...>) {
 *     return (std::string{parts[Is]} + ...);
 *   }
 * }
 * auto c46 = detail_c46::expand(std::make_index_sequence<detail_c46::parts.size()>{});
 *
 */

#include "models/topology/engine.hpp"
#include "models/topology/compute.hpp"

// clang-format off
namespace lahuta::models::topology {
void ModelTopologyEngine::initialize(const ModelTopologyBuildingOptions &opts) {

  auto* atoms_params = get_parameters<ModelAtomsParams>(ModelAtomsComputation<>::label);
  atoms_params->graph_type = opts.graph_type;

  engine_->set_auto_heal(opts.auto_heal);

  engine_->enable(ModelAtomsComputation<>::label,      true);  // required
  engine_->enable(ModelPositionsComputation<>::label,  true);  // required
  engine_->enable(ModelBuildComputation<>::label,      true);  // required

  engine_->enable(ModelBondsComputation<>::label,
                  has_flag(opts.enabled_computations, ModelTopologyComputation::Bonds));
  engine_->enable(ModelAromaticsComputation<>::label,
                  has_flag(opts.enabled_computations, ModelTopologyComputation::Aromatics));
  engine_->enable(ModelDisulfidesComputation<>::label,
                  has_flag(opts.enabled_computations, ModelTopologyComputation::Disulfides));
}

} // namespace lahuta::models::topology
