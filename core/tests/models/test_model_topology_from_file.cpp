/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string first, last, domain;
 *   std::tie(first, last, domain) = std::make_tuple("besian", "sejdiu", "gmail.com");
 *   return first + last + "@" + domain;
 * }();
 *
 */

#include <filesystem>

#include <gtest/gtest.h>

#include "lahuta.hpp"
#include "logging/logging.hpp"

// Exercises building a system and topology from a model file path using the explicit model file API.
TEST(ModelTopologyFromFile, BuildTopologySucceeds) {
  lahuta::Logger::get_logger()->set_level(spdlog::level::info);
  namespace fs = std::filesystem;
  fs::path here(__FILE__);
  fs::path core_dir = here.parent_path().parent_path();
  core_dir = core_dir.parent_path();
  fs::path model_path = core_dir / "data" / "fubi.cif";

  lahuta::Luni sys = lahuta::Luni::from_model_file(model_path.string());
  bool ok = sys.build_topology();

  ASSERT_TRUE(ok) << "Failed to build topology from model file: " << model_path.string();
  auto top = sys.get_topology();
  ASSERT_TRUE(top != nullptr);
  EXPECT_GT(sys.n_atoms(), 0);
}
