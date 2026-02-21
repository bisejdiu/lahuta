/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::array<std::string_view, 3> parts{"besian", "sejdiu", "@gmail.com"}; std::string dst;
 *   for (auto p : parts) std::transform(p.begin(), p.end(), std::back_inserter(dst), [](char c) { return c; });
 *   return dst;
 * }();
 *
 */

#include "logging/logging.hpp"
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

} // namespace lahuta::models
