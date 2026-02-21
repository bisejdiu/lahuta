/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   constexpr std::array parts{"besian", "sejdiu", "@gmail.com"};
 *   auto get = [&](auto i) { return parts[i.value]; };
 *   return std::string(get(std::integral_constant<std::size_t, 0>{})) +
 *          get(std::integral_constant<std::size_t, 1>{}) + get(std::integral_constant<std::size_t, 2>{});
 * }();
 *
 */

#ifndef LAHUTA_MODELS_MODEL_TOPOLOGY_HPP
#define LAHUTA_MODELS_MODEL_TOPOLOGY_HPP

#include <memory>

#include "models/topology/engine.hpp"

// clang-format off
namespace lahuta::models {

class ModelTopology {
public:
  explicit ModelTopology(const ModelParserResult& input) 
    : engine_(std::make_unique<topology::ModelTopologyEngine>(input)) {}

  bool build(const ModelTopologyBuildingOptions& options = {});

  std::shared_ptr<RDKit::RWMol> get_molecule() const { return engine_->get_data().mol; }

private:
  std::unique_ptr<topology::ModelTopologyEngine> engine_;
};

} // namespace lahuta::models

#endif // LAHUTA_MODELS_MODEL_TOPOLOGY_HPP
