// Exercises EnsureTyping end-to-end with the compute engine:
// - Constructs a PipelineContext and runs the following chain using the
//   compute engine with auto-heal enabled:
//     SystemRead -> BuildTopology -> EnsureTyping (Arpeggio or MolStar via dynamic label)
// - It verifies:
//   - The TaskContext sentinel key "atom_typing_mode" is set to the expected
//     provider string ("arpeggio" or "molstar").
//   - The current atom typing mode reported by the topology engine matches
//     the requested mode.
#include <gtest/gtest.h>

#include <filesystem>

#include "analysis/system/computation.hpp"
#include "analysis/topology/computation.hpp"
#include "analysis/topology/ensure_typing.hpp"
#include "compute/engine.hpp"
#include "pipeline/compute/context.hpp"
#include "pipeline/compute/parameters.hpp"

using namespace lahuta;
using namespace lahuta::pipeline::compute;
using lahuta::topology::compute::ComputeEngine;
using lahuta::topology::compute::Mut;

// clang-format off
namespace {

// Run a computation and assert success
static void run_ok(ComputeEngine<PipelineContext, Mut::ReadWrite>& eng, const ComputationLabel& lbl) {
  ASSERT_TRUE(eng.run_from<void>(lbl));
  auto res = eng.get_computation_result(lbl);
  ASSERT_TRUE(res.is_success()) << (res.has_error() ? res.error().get_message() : "unknown error");
}

// Fetch current atom typing mode from topology via engine
static bool current_is_molstar(const Topology& top) {
  auto& eng = const_cast<Topology&>(top).get_engine();
  const auto& lbl = ::lahuta::topology::AtomTypingComputation<>::label;
  auto* p = eng.get_parameters<::lahuta::topology::AtomTypingParams>(lbl);
  return p ? (p->mode == AtomTypingMethod::Molstar) : true;
}

TEST(EnsureTypingTest, SwitchesToArpeggioFromDefaultMolstar) {
  PipelineContext pcx;

  namespace fs = std::filesystem;
  fs::path here(__FILE__);
  fs::path core_dir = here.parent_path().parent_path().parent_path(); // core/tests/.. -> core/
  fs::path item_path = core_dir / "data" / "fubi.cif";

  pcx.item_path = item_path.string();
  lahuta::pipeline::dynamic::TaskContext tctx; pcx.ctx = &tctx;

  // Build engine and register computations
  ComputeEngine<PipelineContext, Mut::ReadWrite> eng(pcx);
  eng.set_auto_heal(true);
  // system
  SystemReadParams sp{}; sp.is_model = false;
  eng.add(std::make_unique<analysis::system::SystemReadComputation>(sp));
  // topology (default typing = Molstar)
  BuildTopologyParams tp{}; tp.flags = TopologyComputation::All; tp.atom_typing_method = AtomTypingMethod::Molstar;
  eng.add(std::make_unique<analysis::topology::BuildTopologyComputation>(tp));
  // ensure arpeggio
  EnsureTypingParams ep{}; ep.desired = AtomTypingMethod::Arpeggio;
  eng.add(std::make_unique<analysis::topology::EnsureTypingComputation>(std::string("ensure_typing_arpeggio"), ep));

  // Execute in order
  run_ok(eng, analysis::system::SystemReadComputation::label);
  run_ok(eng, analysis::topology::BuildTopologyComputation::label);
  run_ok(eng, ComputationLabel{"ensure_typing_arpeggio"});

  auto topo = tctx.get_object<Topology>("topology");
  ASSERT_TRUE(topo);
  EXPECT_FALSE(current_is_molstar(*topo));
  const std::string* s = tctx.get_text("atom_typing_mode");
  ASSERT_NE(s, nullptr);
  EXPECT_EQ(*s, "arpeggio");
}

TEST(EnsureTypingTest, StaysMolstarWhenRequestedMolstar) {
  PipelineContext pcx;
  namespace fs = std::filesystem;
  fs::path here(__FILE__);
  fs::path core_dir = here.parent_path().parent_path().parent_path(); // core/tests/.. -> core/
  fs::path item_path = core_dir / "data" / "1kx2_small.cif";
  pcx.item_path = item_path.string();

  lahuta::pipeline::dynamic::TaskContext tctx; pcx.ctx = &tctx;

  ComputeEngine<PipelineContext, Mut::ReadWrite> eng(pcx);
  eng.set_auto_heal(true);
  SystemReadParams sp{}; sp.is_model = false;
  eng.add(std::make_unique<analysis::system::SystemReadComputation>(sp));
  BuildTopologyParams tp{}; tp.flags = TopologyComputation::All; tp.atom_typing_method = AtomTypingMethod::Molstar;
  eng.add(std::make_unique<analysis::topology::BuildTopologyComputation>(tp));
  EnsureTypingParams ep{}; ep.desired = AtomTypingMethod::Molstar;
  eng.add(std::make_unique<analysis::topology::EnsureTypingComputation>(std::string("ensure_typing_molstar"), ep));

  run_ok(eng, analysis::system::SystemReadComputation::label);
  run_ok(eng, analysis::topology::BuildTopologyComputation::label);
  run_ok(eng, ComputationLabel{"ensure_typing_molstar"});

  auto topo = tctx.get_object<Topology>("topology");
  ASSERT_TRUE(topo);
  EXPECT_TRUE(current_is_molstar(*topo));
  const std::string* s = tctx.get_text("atom_typing_mode");
  ASSERT_NE(s, nullptr);
  EXPECT_EQ(*s, "molstar");
}

} // namespace
