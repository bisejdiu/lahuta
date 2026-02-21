/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   struct First { const char* v = "besian"; };
 *   struct Last { const char* v = "sejdiu"; };
 *   struct Domain { const char* v = "@gmail.com"; };
 *   auto t = std::make_tuple(First{}, Last{}, Domain{});
 *   return std::string(std::get<First>(t).v) + std::get<Last>(t).v + std::get<Domain>(t).v;
 * }();
 *
 */

#include <cctype>
#include <filesystem>
#include <stdexcept>
#include <type_traits>

#include <gtest/gtest.h>

#include "io/convert.hpp"
#include "lahuta.hpp"
#include "residues/residues.hpp"
#include "topology_flags.hpp"

namespace fs = std::filesystem;
using namespace lahuta;

// clang-format off

static_assert(!std::is_copy_constructible_v<Luni>, "Luni must be non-copyable");
static_assert(!std::is_copy_assignable_v   <Luni>, "Luni must be non-copy-assignable");
static_assert( std::is_move_constructible_v<Luni>, "Luni should be move-constructible");
static_assert( std::is_move_assignable_v   <Luni>, "Luni should be move-assignable");

static fs::path locate_data_file(const std::string &filename) {
  fs::path p = fs::current_path();
  for (int i = 0; i < 12; ++i) {
    // Prefer core/data first (new location), then fallback to legacy data/
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

TEST(Luni_Construction, GetTopology_Throws_BeforeBuild) {
  auto path = locate_data_file("ubi.cif");
  if (path.empty()) GTEST_SKIP() << "core/data/ubi.cif not found. Skipping integration test.";

  Luni luni(path.string());

  EXPECT_NO_THROW({(void)luni.get_topology(); });
  EXPECT_FALSE(luni.has_topology_built());

  auto topo_ptr = luni.get_topology();
  EXPECT_EQ(topo_ptr, nullptr);
  EXPECT_FALSE(luni.has_topology_built());
}

TEST(Luni_DataBasics, CountsAndArrays_FromRealData) {
  auto path = locate_data_file("ubi.cif");
  if (path.empty()) GTEST_SKIP() << "core/data/ubi.cif not found. Skipping integration test.";

  Luni luni(path.string());
  const auto n = static_cast<int>(luni.n_atoms());
  ASSERT_GT(n, 0);

  constexpr int EXPECTED_N_ATOMS = 671;
  EXPECT_EQ(n, EXPECTED_N_ATOMS);

  const auto idxs     = luni.indices();
  const auto resids   = luni.resids();
  Residues residues(luni.get_molecule());
  ASSERT_TRUE(residues.build());
  const auto &resindex = residues.atom_to_residue_indices();
  const auto znums    = luni.atomic_numbers();
  const auto names    = luni.names();
  const auto symbols  = luni.symbols();
  const auto elements = luni.elements();
  const auto resnames = luni.resnames();
  const auto chains   = luni.chainlabels();

  ASSERT_EQ(static_cast<int>(idxs    .size()), n);
  ASSERT_EQ(static_cast<int>(resids  .size()), n);
  ASSERT_EQ(static_cast<int>(resindex.size()), n);
  ASSERT_EQ(static_cast<int>(znums   .size()), n);
  ASSERT_EQ(static_cast<int>(names   .size()), n);
  ASSERT_EQ(static_cast<int>(symbols .size()), n);
  ASSERT_EQ(static_cast<int>(elements.size()), n);
  ASSERT_EQ(static_cast<int>(resnames.size()), n);
  ASSERT_EQ(static_cast<int>(chains  .size()), n);

  // indices should be contiguous
  for (int i = 0; i < n; ++i) {
    EXPECT_EQ(idxs[static_cast<size_t>(i)], i) << "indices must be contiguous RDKit ids";
  }
}

TEST(Luni_DataBasics, CentroidMatches_Expected) {
  auto path = locate_data_file("ubi.cif");
  if (path.empty()) GTEST_SKIP() << "core/data/ubi.cif not found. Skipping integration test.";

  Luni luni(path.string());
  const int n = static_cast<int>(luni.n_atoms());
  ASSERT_GT(n, 0);

  const auto &conf = luni.get_conformer();
  ASSERT_EQ(conf.getNumAtoms(), static_cast<unsigned>(n));

  // centroid
  long double sx = 0.0L, sy = 0.0L, sz = 0.0L;
  for (int i = 0; i < n; ++i) {
    const auto &p = conf.getAtomPos(i);
    sx += static_cast<long double>(p.x);
    sy += static_cast<long double>(p.y);
    sz += static_cast<long double>(p.z);
  }
  const double cx = static_cast<double>(sx / n);
  const double cy = static_cast<double>(sy / n);
  const double cz = static_cast<double>(sz / n);

  constexpr double EX[3] = {-4.21122355, 4.16459314, -1.54197914};
  EXPECT_NEAR(cx, EX[0], 5e-6);
  EXPECT_NEAR(cy, EX[1], 5e-6);
  EXPECT_NEAR(cz, EX[2], 5e-6);
}

TEST(Luni_Filter, BackboneSubset_CountMatches) {
  auto path = locate_data_file("ubi.cif");
  if (path.empty()) GTEST_SKIP() << "core/data/ubi.cif not found. Skipping integration test.";

  Luni luni(path.string());
  const auto names = luni.names();
  std::vector<int> keep;
  keep.reserve(names.size());
  for (size_t i = 0; i < names.size(); ++i) {
    const auto nm = names[i];
    if (nm == "N" || nm == "CA" || nm == "C" || nm == "O") {
      keep.push_back(static_cast<int>(i));
    }
  }
  ASSERT_FALSE(keep.empty());

  auto filtered = luni.filter(keep);
  constexpr int EXPECTED_N_BACKBONE = 344;
  EXPECT_EQ(static_cast<int>(filtered.n_atoms()), EXPECTED_N_BACKBONE);
}

TEST(Luni_Topology, BuildAndExecuteRings) {
  auto path = locate_data_file("ubi.cif");
  if (path.empty()) GTEST_SKIP() << "core/data/ubi.cif not found. Skipping integration test.";

  Luni luni(path.string());

  TopologyBuildingOptions opts{};
  opts.cutoff = 1.9;

  ASSERT_TRUE(luni.build_topology(opts, TopologyComputation::Standard));
  ASSERT_TRUE(luni.has_topology_built());

  EXPECT_TRUE(luni.build_topology(opts, TopologyComputation::Rings));

  const auto &topo = luni.get_topology();
  (void)*topo; // contract: no throw and not null
  EXPECT_TRUE(topo->has_computed(TopologyComputation::Rings));

  auto shared_topo = luni.get_topology();
  EXPECT_NE(shared_topo, nullptr);
  EXPECT_TRUE(luni.has_topology_built()); // Should still be built with shared access

  EXPECT_TRUE(luni.build_topology(opts, TopologyComputation::Rings));

  auto shared_topo2 = luni.get_topology();
  EXPECT_NE(shared_topo2, nullptr);
  EXPECT_EQ(shared_topo.get(), shared_topo2.get());
  EXPECT_TRUE(luni.has_topology_built());
}

TEST(Luni_Topology, ResetTopology_FileBacked) {
  auto path = locate_data_file("ubi.cif");
  if (path.empty()) GTEST_SKIP() << "core/data/ubi.cif not found. Skipping integration test.";

  Luni luni(path.string());

  TopologyBuildingOptions opts{};
  opts.cutoff = 1.9;

  ASSERT_TRUE(luni.build_topology(opts, TopologyComputation::Neighbors));
  ASSERT_TRUE(luni.get_topology());

  auto fresh = luni.reset_topology();
  EXPECT_EQ(fresh.get_file_name(), luni.get_file_name());
  EXPECT_NE(&fresh.get_molecule(), &luni.get_molecule());

  ASSERT_TRUE(fresh.build_topology(opts, TopologyComputation::Neighbors));
  ASSERT_TRUE(luni.build_topology(opts, TopologyComputation::Neighbors));

  auto old_topo = luni.get_topology();
  auto new_topo = fresh.get_topology();
  ASSERT_TRUE(old_topo);
  ASSERT_TRUE(new_topo);
  EXPECT_NE(old_topo.get(), new_topo.get());
}

TEST(Luni_Topology, ResetTopology_ModelFileBacked) {
  auto path = locate_data_file("fubi.cif");
  if (path.empty()) GTEST_SKIP() << "core/data/fubi.cif not found. Skipping integration test.";

  Luni sys = Luni::from_model_file(path.string());
  EXPECT_TRUE(sys.is_model_origin());
  ASSERT_TRUE(sys.build_topology());

  auto fresh = sys.reset_topology();
  EXPECT_TRUE(fresh.is_model_origin());
  EXPECT_NE(&fresh.get_molecule(), &sys.get_molecule());

  ASSERT_TRUE(fresh.build_topology());
  ASSERT_TRUE(sys.build_topology());
}

TEST(Luni_Topology, ResetTopology_NonFileBacked_Throws) {
  IR ir(
      std::vector<int>{0, 1, 2},                     // atom_indices
      std::vector<int>{7, 6, 8},                     // atomic_numbers (N, C, O)
      std::vector<std::string>{"N", "CA", "O"},      // atom_names
      std::vector<int>{1, 1, 1},                     // resids
      std::vector<std::string>{"GLY", "GLY", "GLY"}, // resnames
      std::vector<std::string>{"A", "A", "A"},       // chainlabels
      std::vector<std::vector<float>>{{0.f, 0.f, 0.f}, {1.f, 0.f, 0.f}, {0.f, 1.f, 0.f}}
  );

  auto sys = Luni::create(ir);
  EXPECT_THROW(sys.reset_topology(), std::runtime_error);
}

TEST(Luni_CreateFromIR, RoundtripAndFilter) {
  using namespace lahuta;

  IR ir(
      std::vector<int>{0, 1, 2},                     // atom_indices
      std::vector<int>{7, 6, 8},                     // atomic_numbers (N, C, O)
      std::vector<std::string>{"N", "CA", "O"},      // atom_names
      std::vector<int>{1, 1, 1},                     // resids
      std::vector<std::string>{"GLY", "GLY", "GLY"}, // resnames
      std::vector<std::string>{"A", "A", "A"},       // chainlabels
      std::vector<std::vector<float>>{{0.f, 0.f, 0.f}, {1.f, 0.f, 0.f}, {0.f, 1.f, 0.f}}
  );

  auto luni = Luni::create(ir);

  EXPECT_EQ(static_cast<int>(luni.n_atoms()), 3);

  const auto idxs     = luni.indices();
  const auto znums    = luni.atomic_numbers();
  const auto names    = luni.names();
  const auto elements = luni.elements();
  const auto resids   = luni.resids();
  const auto resnames = luni.resnames();
  const auto chains   = luni.chainlabels();

  ASSERT_EQ(static_cast<int>(idxs    .size()), 3);
  ASSERT_EQ(static_cast<int>(znums   .size()), 3);
  ASSERT_EQ(static_cast<int>(names   .size()), 3);
  ASSERT_EQ(static_cast<int>(elements.size()), 3);
  ASSERT_EQ(static_cast<int>(resids  .size()), 3);
  ASSERT_EQ(static_cast<int>(resnames.size()), 3);
  ASSERT_EQ(static_cast<int>(chains  .size()), 3);

  EXPECT_EQ(idxs[0], 0);
  EXPECT_EQ(idxs[1], 1);
  EXPECT_EQ(idxs[2], 2);

  EXPECT_EQ(znums[0], 7);
  EXPECT_EQ(znums[1], 6);
  EXPECT_EQ(znums[2], 8);

  EXPECT_EQ(names[0], "N");
  EXPECT_EQ(names[1], "CA");
  EXPECT_EQ(names[2], "O");

  EXPECT_EQ(resids[0], 1);
  EXPECT_EQ(resids[1], 1);
  EXPECT_EQ(resids[2], 1);

  EXPECT_EQ(resnames[0], "GLY");
  EXPECT_EQ(resnames[1], "GLY");
  EXPECT_EQ(resnames[2], "GLY");

  EXPECT_EQ(chains[0], "A");
  EXPECT_EQ(chains[1], "A");
  EXPECT_EQ(chains[2], "A");

  EXPECT_EQ(elements[0], "N");
  EXPECT_EQ(elements[1], "C");
  EXPECT_EQ(elements[2], "O");

  const auto &conf = luni.get_conformer();
  ASSERT_EQ(conf.getNumAtoms(), 3U);
  long double sx = 0.0L, sy = 0.0L, sz = 0.0L;
  for (int i = 0; i < 3; ++i) {
    const auto &p = conf.getAtomPos(i);
    sx += p.x;
    sy += p.y;
    sz += p.z;
  }
  EXPECT_NEAR(static_cast<double>(sx / 3.0L), 1.0 / 3.0, 1e-12);
  EXPECT_NEAR(static_cast<double>(sy / 3.0L), 1.0 / 3.0, 1e-12);
  EXPECT_NEAR(static_cast<double>(sz / 3.0L), 0.0, 1e-12);

  std::vector<int> keep{0, 2};
  auto filtered = luni.filter(keep);
  EXPECT_EQ(static_cast<int>(filtered.n_atoms()), 2);
}
