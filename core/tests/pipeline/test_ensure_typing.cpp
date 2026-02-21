/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s;
 *   return (s += "besian", s += "sejdiu", s += "@gmail.com", s);
 * }();
 *
 */

// Exercises EnsureTyping end-to-end with the compute engine:
// - Constructs a PipelineContext and runs the following chain using the
//   compute engine with auto-heal enabled:
//     SystemRead -> BuildTopology -> EnsureTyping (Arpeggio or MolStar via dynamic label)
// - It verifies:
//   - The TaskContext sentinel key "atom_typing_mode" is set to the expected
//     provider string ("arpeggio" or "molstar").
//   - The current atom typing mode reported by the topology engine matches
//     the requested mode.
#include <filesystem>

#include <gtest/gtest.h>

#include "analysis/system/computation.hpp"
#include "analysis/topology/computation.hpp"
#include "analysis/topology/ensure_typing.hpp"
#include "compute/engine.hpp"
#include "pipeline/task/compute/context.hpp"
#include "pipeline/task/compute/parameters.hpp"

using namespace lahuta;
namespace P = lahuta::pipeline;
namespace C = lahuta::compute;

namespace {

#if defined(__has_feature)
#  if __has_feature(address_sanitizer)
constexpr bool RunningUnderASan = true;
#  else
constexpr bool RunningUnderASan = false;
#  endif
#elif defined(__SANITIZE_ADDRESS__)
constexpr bool RunningUnderASan = true;
#else
constexpr bool RunningUnderASan = false;
#endif

// Run a computation and assert success
static void run_ok(C::ComputeEngine<P::PipelineContext> &eng, const C::ComputationLabel &lbl) {
  ASSERT_TRUE(eng.run_from<void>(lbl));
  auto res = eng.get_computation_result(lbl);
  ASSERT_TRUE(res.is_success()) << (res.has_error() ? res.error().get_message() : "unknown error");
}

// Fetch current atom typing mode from topology via engine
static bool current_is_molstar(const Topology &top) {
  auto &eng       = const_cast<Topology &>(top).get_engine();
  const auto &lbl = ::lahuta::topology::AtomTypingComputation<>::label;
  auto *p         = eng.get_parameters<::lahuta::topology::AtomTypingParams>(lbl);
  return p ? (p->mode == AtomTypingMethod::Molstar) : true;
}

TEST(EnsureTypingTest, SwitchesToArpeggioFromDefaultMolstar) {
  P::PipelineContext pcx;

  namespace fs = std::filesystem;
  fs::path here(__FILE__);
  fs::path core_dir  = here.parent_path().parent_path().parent_path(); // core/tests/.. -> core/
  fs::path item_path = core_dir / "data" / "fubi.cif";

  pcx.item_path = item_path.string();
  P::TaskContext tctx;
  pcx.ctx = &tctx;

  // Build engine and register computations
  C::ComputeEngine<P::PipelineContext> eng(pcx);
  // system
  P::SystemReadParams sp{};
  sp.is_model = false;
  eng.add(std::make_unique<analysis::SystemReadComputation>(sp));
  // topology (default typing = Molstar)
  P::BuildTopologyParams tp{};
  tp.flags              = TopologyComputation::All;
  tp.atom_typing_method = AtomTypingMethod::Molstar;
  eng.add(std::make_unique<analysis::BuildTopologyComputation>(tp));
  // ensure arpeggio
  P::EnsureTypingParams ep{};
  ep.desired = AtomTypingMethod::Arpeggio;
  eng.add(std::make_unique<analysis::EnsureTypingComputation>(std::string("ensure_typing_arpeggio"), ep));

  // Execute in order
  run_ok(eng, analysis::SystemReadComputation::label);
  run_ok(eng, analysis::BuildTopologyComputation::label);
  run_ok(eng, C::ComputationLabel{"ensure_typing_arpeggio"});

  auto topo = tctx.topology();
  ASSERT_TRUE(topo);
  EXPECT_FALSE(current_is_molstar(*topo));
  const std::string *s = tctx.get_text("atom_typing_mode");
  ASSERT_NE(s, nullptr);
  EXPECT_EQ(*s, "Arpeggio");
}

TEST(EnsureTypingTest, StaysMolstarWhenRequestedMolstar) {
  if (RunningUnderASan) {
    GTEST_SKIP() << "Ring perception by RDKit is flaky when AddressSanitizer is enabled.";
  }

  P::PipelineContext pcx;
  namespace fs = std::filesystem;
  fs::path here(__FILE__);
  fs::path core_dir  = here.parent_path().parent_path().parent_path(); // core/tests/.. -> core/
  fs::path item_path = core_dir / "data" / "1kx2_small.cif";
  pcx.item_path      = item_path.string();

  P::TaskContext tctx;
  pcx.ctx = &tctx;

  C::ComputeEngine<P::PipelineContext> eng(pcx);
  P::SystemReadParams sp{};
  sp.is_model = false;
  eng.add(std::make_unique<analysis::SystemReadComputation>(sp));
  P::BuildTopologyParams tp{};
  tp.flags              = TopologyComputation::All;
  tp.atom_typing_method = AtomTypingMethod::Molstar;
  eng.add(std::make_unique<analysis::BuildTopologyComputation>(tp));
  P::EnsureTypingParams ep{};
  ep.desired = AtomTypingMethod::Molstar;
  eng.add(std::make_unique<analysis::EnsureTypingComputation>(std::string("ensure_typing_molstar"), ep));

  run_ok(eng, analysis::SystemReadComputation::label);
  run_ok(eng, analysis::BuildTopologyComputation::label);
  run_ok(eng, C::ComputationLabel{"ensure_typing_molstar"});

  auto topo = tctx.topology();
  ASSERT_TRUE(topo);
  EXPECT_TRUE(current_is_molstar(*topo));
  const std::string *s = tctx.get_text("atom_typing_mode");
  ASSERT_NE(s, nullptr);
  EXPECT_EQ(*s, "MolStar");
}

} // namespace
