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

#include <array>
#include <filesystem>
#include <random>

#include <gtest/gtest.h>

#include "lahuta.hpp"
#include "topology/compute.hpp"
#include "topology_flags.hpp"

namespace fs = std::filesystem;
using namespace lahuta;

namespace {

static fs::path locate_data_file(const std::string &filename) {
  fs::path p = fs::current_path();
  for (int i = 0; i < 12; ++i) {
    fs::path cand_core = p / "core" / "data" / filename;
    if (fs::exists(cand_core)) return cand_core;

    fs::path cand_legacy = p / "data" / filename;
    if (fs::exists(cand_legacy)) return cand_legacy;

    if (p.has_parent_path()) {
      p = p.parent_path();
    } else {
      break;
    }
  }
  return {};
}

struct CompInfo {
  TopologyComputation flag;
  const topology::ComputationLabel &label;
};

static constexpr std::array<CompInfo, 6> Base = {
    {
     {TopologyComputation::Neighbors, topology::NeighborSearchComputation<>::label},
     {TopologyComputation::Bonds, topology::BondComputation<>::label},
     {TopologyComputation::NonStandardBonds, topology::NonStandardBondComputation<>::label},
     {TopologyComputation::Residues, topology::ResidueComputation<>::label},
     {TopologyComputation::Rings, topology::RingComputation<>::label},
     {TopologyComputation::AtomTyping, topology::AtomTypingComputation<>::label},
     }
};

static TopologyComputation expand_dependencies(TopologyComputation mask) {
  bool changed = true;
  while (changed) {
    changed = false;

    if (has_flag(mask, TopologyComputation::AtomTyping) && !has_flag(mask, TopologyComputation::Rings)) {
      mask    = mask | TopologyComputation::Rings;
      changed = true;
    }
    if (has_flag(mask, TopologyComputation::Rings)) {
      if (!has_flag(mask, TopologyComputation::Bonds)) {
        mask    = mask | TopologyComputation::Bonds;
        changed = true;
      }
      if (!has_flag(mask, TopologyComputation::Residues)) {
        mask    = mask | TopologyComputation::Residues;
        changed = true;
      }
    }
    if (has_flag(mask, TopologyComputation::NonStandardBonds) &&
        !has_flag(mask, TopologyComputation::Bonds)) {
      mask    = mask | TopologyComputation::Bonds;
      changed = true;
    }
    if (has_flag(mask, TopologyComputation::Bonds) && !has_flag(mask, TopologyComputation::Neighbors)) {
      mask    = mask | TopologyComputation::Neighbors;
      changed = true;
    }
  }
  return mask;
}

} // namespace

TEST(TopologyIdempotency, RepeatedRequestsDoNotRecompute) {
  auto path = locate_data_file("1kx2_small.cif");
  if (path.empty()) GTEST_SKIP() << "core/data/1kx2_small.cif not found. Skipping test.";

  Luni luni(path.string());
  TopologyBuildingOptions opts{};
  opts.compute_nonstandard_bonds = true;
  opts.atom_typing_method        = AtomTypingMethod::Molstar;
  opts.mode                      = TopologyBuildMode::Generic;

  const std::array<TopologyComputation, 10> candidates = {
      {
       TopologyComputation::None,
       TopologyComputation::Neighbors,
       TopologyComputation::Bonds,
       TopologyComputation::NonStandardBonds,
       TopologyComputation::Residues,
       TopologyComputation::Rings,
       TopologyComputation::AtomTyping,
       TopologyComputation::Standard,
       TopologyComputation::Extended,
       TopologyComputation::Complete,
       }
  };

  std::mt19937 rng(1337);
  std::uniform_int_distribution<std::size_t> dist(0, candidates.size() - 1);

  TopologyComputation ever_required = TopologyComputation::None;

  for (int i = 0; i < 64; ++i) {
    auto include  = candidates[dist(rng)];
    ever_required = ever_required | expand_dependencies(include);

    ASSERT_TRUE(luni.build_topology(opts, include));
    auto topo = luni.get_topology();
    ASSERT_TRUE(topo);
    auto *eng = topo->get_engine().get_engine();
    ASSERT_NE(eng, nullptr);

    for (const auto &info : Base) {
      const auto count = eng->get_run_count(info.label);
      EXPECT_LE(count, 1u);
      if (has_flag(ever_required, info.flag)) {
        EXPECT_GE(count, 1u);
      } else {
        EXPECT_EQ(count, 0u);
      }
    }

    EXPECT_TRUE(topo->has_computed(ever_required));
  }
}

TEST(TopologyIdempotency, NoRecomputeAfterFullBuild) {
  auto path = locate_data_file("1kx2_small.cif");
  if (path.empty()) GTEST_SKIP() << "core/data/1kx2_small.cif not found. Skipping test.";

  Luni luni(path.string());
  TopologyBuildingOptions opts{};
  opts.compute_nonstandard_bonds = true;
  opts.atom_typing_method        = AtomTypingMethod::Molstar;
  opts.mode                      = TopologyBuildMode::Generic;

  ASSERT_TRUE(luni.build_topology(opts, TopologyComputation::Complete));
  auto topo = luni.get_topology();
  ASSERT_TRUE(topo);
  auto *eng = topo->get_engine().get_engine();
  ASSERT_NE(eng, nullptr);

  eng->reset_run_counts();

  const std::array<TopologyComputation, 10> candidates = {
      {
       TopologyComputation::None,
       TopologyComputation::Neighbors,
       TopologyComputation::Bonds,
       TopologyComputation::NonStandardBonds,
       TopologyComputation::Residues,
       TopologyComputation::Rings,
       TopologyComputation::AtomTyping,
       TopologyComputation::Standard,
       TopologyComputation::Extended,
       TopologyComputation::Complete,
       }
  };

  std::mt19937 rng(4242);
  std::uniform_int_distribution<std::size_t> dist(0, candidates.size() - 1);

  for (int i = 0; i < 128; ++i) {
    auto include = candidates[dist(rng)];
    ASSERT_TRUE(luni.build_topology(opts, include));

    for (const auto &info : Base) {
      EXPECT_EQ(eng->get_run_count(info.label), 0u);
    }
  }
}

TEST(TopologyIdempotency, NoneDoesNotExecute) {
  auto path = locate_data_file("1kx2_small.cif");
  if (path.empty()) GTEST_SKIP() << "core/data/1kx2_small.cif not found. Skipping test.";

  Luni luni(path.string());
  TopologyBuildingOptions opts{};
  opts.compute_nonstandard_bonds = true;
  opts.atom_typing_method        = AtomTypingMethod::Molstar;
  opts.mode                      = TopologyBuildMode::Generic;

  ASSERT_TRUE(luni.build_topology(opts, TopologyComputation::None));
  auto topo = luni.get_topology();
  ASSERT_TRUE(topo);
  auto *eng = topo->get_engine().get_engine();
  ASSERT_NE(eng, nullptr);

  for (const auto &info : Base) {
    EXPECT_EQ(eng->get_run_count(info.label), 0u);
    EXPECT_FALSE(topo->has_computed(info.flag));
  }

  for (int i = 0; i < 8; ++i) {
    ASSERT_TRUE(luni.build_topology(opts, TopologyComputation::None));
  }

  for (const auto &info : Base) {
    EXPECT_EQ(eng->get_run_count(info.label), 0u);
    EXPECT_FALSE(topo->has_computed(info.flag));
  }
}
