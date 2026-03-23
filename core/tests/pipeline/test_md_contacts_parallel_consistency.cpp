/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s;
 *   s += "besian";
 *   s += "sejdiu";
 *   s += "@gmail.com";
 *   return s;
 * }();
 *
 */

#include <algorithm>
#include <filesystem>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "analysis/contacts/computation.hpp"
#include "pipeline/ingest/descriptor.hpp"
#include "pipeline/runtime/api.hpp"
#include "pipeline/task/api.hpp"
#include "sinks/memory.hpp"

using namespace lahuta;
namespace P = lahuta::pipeline;

namespace {

namespace fs = std::filesystem;

static fs::path locate_simulation_file(const std::string &filename) {
#ifdef LAHUTA_SIMDB_DIR
  {
    fs::path base{LAHUTA_SIMDB_DIR};
    fs::path cand = base / filename;
    if (fs::exists(cand)) return cand;
  }
#endif
  fs::path p = fs::current_path();
  for (int i = 0; i < 12; ++i) {
    fs::path cand_core = p / "core" / "data" / "simulationdatabase" / filename;
    if (fs::exists(cand_core)) return cand_core;

    if (!p.has_parent_path()) break;
    p = p.parent_path();
  }
  return {};
}

struct SingleTrajectoryDescriptor final : P::IDescriptor {
  std::string structure;
  std::vector<std::string> xtcs;
  bool done{false};

  SingleTrajectoryDescriptor(std::string s, std::vector<std::string> x)
      : structure(std::move(s)), xtcs(std::move(x)) {}

  std::optional<P::IngestDescriptor> next() override {
    if (done) return std::nullopt;
    done = true;
    P::IngestDescriptor d;
    d.id     = structure;
    d.origin = P::MDRef{structure, xtcs};
    return d;
  }

  void reset() override { done = false; }
};

struct ContactsRun {
  P::StageManager::RunReport report;
  std::vector<std::string> payloads;
};

static ContactsRun run_md_contacts(const fs::path &structure_path, const fs::path &trajectory_path,
                                   std::size_t trajectory_repeats, std::size_t threads) {
  std::vector<std::string> xtcs(trajectory_repeats, trajectory_path.string());
  auto src = std::make_shared<SingleTrajectoryDescriptor>(structure_path.string(), std::move(xtcs));

  P::StageManager mgr(std::move(src));
  mgr.set_auto_builtins(true);

  P::ContactsParams p{};
  p.provider = analysis::ContactProvider::MolStar;
  p.channel  = "contacts";
  p.format   = P::ContactsOutputFormat::Json;

  mgr.add_computation(
      "contacts",
      {},
      [label = std::string("contacts"), p]() {
        return std::make_unique<analysis::ContactsComputation>(label, p);
      },
      /*thread_safe=*/true);

  auto sink = std::make_shared<P::MemorySink>();
  mgr.connect_sink("contacts", sink);

  mgr.compile();
  auto report   = mgr.run(threads);
  auto payloads = sink->result();
  std::sort(payloads.begin(), payloads.end());

  return ContactsRun{std::move(report), std::move(payloads)};
}

} // namespace

TEST(MdContactsParallelConsistency, ParallelOutputMatchesSingleThreadedBaseline) {
  const auto structure_path  = locate_simulation_file("msd_coords.gro");
  const auto trajectory_path = locate_simulation_file("msd_traj.xtc");
  if (structure_path.empty() || trajectory_path.empty()) {
    GTEST_SKIP() << "Simulation data not available";
  }

  // Repeat the tiny fixture enough times to create meaningful concurrency while
  // keeping the test fast.
  constexpr std::size_t kTrajectoryRepeats = 32;

  const auto baseline = run_md_contacts(structure_path, trajectory_path, kTrajectoryRepeats, 1);
  const auto parallel = run_md_contacts(structure_path, trajectory_path, kTrajectoryRepeats, 8);

  ASSERT_EQ(baseline.report.items_skipped, std::size_t{0});
  ASSERT_EQ(parallel.report.items_skipped, std::size_t{0});
  ASSERT_EQ(baseline.report.items_processed, parallel.report.items_processed);
  ASSERT_EQ(baseline.payloads.size(), baseline.report.items_processed);
  ASSERT_EQ(parallel.payloads.size(), parallel.report.items_processed);
  EXPECT_EQ(parallel.payloads, baseline.payloads);
}
