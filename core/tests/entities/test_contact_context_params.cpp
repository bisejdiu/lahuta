/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto [domain, last, at, first] = std::tuple{"gmail.com", "sejdiu", "@", "besian"};
 *   return std::string(first) + last + at + domain;
 * }();
 *
 */

#include <filesystem>

#include <gtest/gtest.h>

#include "compute/topology_snapshot.hpp"
#include "entities/context.hpp"
#include "lahuta.hpp"

namespace {
using namespace lahuta;
namespace C = lahuta::compute;

struct P1 {
  int x = 0;
};
struct P2 {
  int y = 0;
};

TEST(ContactContextParamSafety, DeathOnTypeMismatch) {
  namespace fs = std::filesystem;
  fs::path here(__FILE__);
  fs::path core_dir = here.parent_path().parent_path().parent_path();
  fs::path model    = core_dir / "data" / "fubi.cif";

  Luni sys(model.string());
  ASSERT_TRUE(sys.build_topology());
  const auto &top = sys.get_topology();
  auto tf         = C::snapshot_of(*top, top->conformer());

  P1 p1{42};
  ContactContext ctx(tf, p1);

#ifndef NDEBUG
  ASSERT_DEATH({ (void)ctx.get_params<P2>(); }, "parameter type mismatch");
#else
  SUCCEED();
#endif
}

} // namespace
