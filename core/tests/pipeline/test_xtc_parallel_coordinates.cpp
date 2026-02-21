/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   constexpr auto f = []() constexpr { return "besian"; };
 *   constexpr auto l = []() constexpr { return "sejdiu"; };
 *   constexpr auto d = []() constexpr { return "@gmail.com"; };
 *   return std::string(f()) + l() + d();
 * }();
 *
 */

#include <atomic>
#include <filesystem>
#include <mutex>
#include <optional>
#include <random>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <GraphMol/Conformer.h>
#include <gtest/gtest.h>

#include "analysis/system/computation.hpp"
#include "compute/dependency.hpp"
#include "compute/result.hpp"
#include "logging/logging.hpp"
#include "pipeline/data/ingestion.hpp"
#include "pipeline/ingest/descriptor.hpp"
#include "pipeline/runtime/manager.hpp"
#include "pipeline/task/compute/context.hpp"

using namespace lahuta;
namespace C = lahuta::compute;
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

static fs::path pdb_path, xtc_path;

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

static double sq_distance(const RDGeom::Point3D &a, const RDGeom::Point3D &b) {
  double dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
  return dx * dx + dy * dy + dz * dz;
}

struct CoordinateCollector {
  std::mutex mtx;
  std::once_flag indices_once;
  std::vector<unsigned> selected_indices; // size <= 100, initialized once
  std::unordered_map<std::size_t, std::vector<RDGeom::Point3D>> coords_by_conformer_id;
  std::atomic<size_t> frames{0};

  void reset() {
    std::scoped_lock lk(mtx);
    selected_indices.clear();
    coords_by_conformer_id.clear();
    frames = 0;
  }
};

static std::vector<unsigned> select_random_unique_indices(std::size_t num_atoms, std::size_t k,
                                                          uint64_t seed) {
  std::vector<unsigned> idxs;
  if (num_atoms == 0 || k == 0) return idxs;
  const std::size_t pick = std::min<std::size_t>(k, num_atoms);
  // Use a deterministic PRNG with fixed seed
  std::mt19937_64 rng(seed);
  std::uniform_int_distribution<std::size_t> dist(0, num_atoms - 1);
  std::unordered_set<std::size_t> used;
  used.reserve(pick * 2);
  while (idxs.size() < pick) {
    std::size_t v = dist(rng);
    if (used.insert(v).second) idxs.push_back(static_cast<unsigned>(v));
  }
  return idxs;
}

struct CoordsCaptureParams : P::ParameterBase<CoordsCaptureParams> {
  static constexpr P::ParameterInterface::TypeId TYPE_ID = 204; // local test id
};

class CoordsCaptureComputation final
    : public P::ReadWriteComputation<P::PipelineContext, CoordsCaptureParams, CoordsCaptureComputation> {
public:
  using Base = P::ReadWriteComputation<P::PipelineContext, CoordsCaptureParams, CoordsCaptureComputation>;
  explicit CoordsCaptureComputation(std::shared_ptr<CoordinateCollector> collector)
      : Base(CoordsCaptureParams{}), collector_(std::move(collector)) {}

  static constexpr P::ComputationLabel label{"coords_capture"};
  using dependencies = C::Dependencies<C::Dependency<analysis::SystemReadComputation, void>>;

  P::ComputationResult execute_typed(P::DataContext<P::PipelineContext> &ctx, const CoordsCaptureParams &) {
    auto &data = ctx.data();

    if (!data.frame) return P::ComputationResult(true);

    // Build a conformer view bound to this frame's coordinates
    auto view = data.frame->load_coordinates();
    auto slab = view.shared_positions ? std::const_pointer_cast<RDGeom::POINT3D_VECT>(view.shared_positions)
                                      : std::make_shared<RDGeom::POINT3D_VECT>(std::move(view.positions));
    RDKit::Conformer conf;
    conf.set3D(true);
    conf.bindExternalPositions(std::move(slab));
    const RDKit::Conformer &cref = conf; // const-view to avoid mutating accessors

    std::call_once(collector_->indices_once, [&] {
      const std::size_t natoms     = static_cast<std::size_t>(cref.getNumAtoms());
      collector_->selected_indices = select_random_unique_indices(natoms, 100, /*seed=*/8351);
    });

    std::vector<RDGeom::Point3D> coords;
    coords.reserve(collector_->selected_indices.size());
    for (unsigned idx : collector_->selected_indices) {
      coords.push_back(cref.getAtomPos(idx));
    }
    {
      std::scoped_lock lk(collector_->mtx);
      collector_->coords_by_conformer_id.emplace(data.conformer_id, std::move(coords));
    }
    ++collector_->frames;
    return P::ComputationResult(true);
  }

private:
  std::shared_ptr<CoordinateCollector> collector_;
};

} // namespace

TEST(ParallelXtcCoordinates, ConsistencyAcrossThreads) {
  Logger::get_instance().set_log_level(Logger::LogLevel::Info);

  pdb_path = locate_simulation_file("lysozyme.gro");
  xtc_path = locate_simulation_file("lysozyme.xtc");
  if (pdb_path.empty() || xtc_path.empty()) {
    GTEST_SKIP() << "Simulation data not available";
  }

  // Oracle
  auto baseline = std::make_shared<CoordinateCollector>();
  baseline->reset();
  {
    auto src = std::make_unique<SingleTrajectoryDescriptor>(std::string(pdb_path),
                                                            std::vector<std::string>{std::string(xtc_path)});
    P::StageManager mgr(std::move(src));
    mgr.set_auto_builtins(true);
    mgr.add_computation("coords_capture", {"system"}, [baseline] {
      return std::make_unique<CoordsCaptureComputation>(baseline);
    });
    mgr.run(1);
  }

  ASSERT_GT(baseline->frames.load(), 0u) << "No frames processed in baseline";
  ASSERT_FALSE(baseline->selected_indices.empty()) << "No indices selected in baseline";

  const auto baseline_indices = baseline->selected_indices; // copy for equality checks
  const auto baseline_frames  = baseline->frames.load();
  const auto &baseline_map    = baseline->coords_by_conformer_id;

  const std::vector<int> threads = {4, 8, 16};
  for (int t : threads) {
    auto trial = std::make_shared<CoordinateCollector>();
    trial->reset();

    auto src = std::make_unique<SingleTrajectoryDescriptor>(std::string(pdb_path),
                                                            std::vector<std::string>{std::string(xtc_path)});
    P::StageManager mgr(std::move(src));
    mgr.set_auto_builtins(true);
    mgr.add_computation("coords_capture", {"system"}, [trial] {
      return std::make_unique<CoordsCaptureComputation>(trial);
    });
    mgr.run(static_cast<std::size_t>(t));

    // Basic invariants
    ASSERT_EQ(trial->frames.load(), baseline_frames)
        << "Thread count " << t << " produced different frame count";
    ASSERT_EQ(trial->selected_indices, baseline_indices) << "Selected indices differ across runs";
    ASSERT_EQ(trial->coords_by_conformer_id.size(), baseline_map.size()) << "Map size mismatch for t=" << t;

    const double tol = 1e-12;
    for (const auto &kv : baseline_map) {
      auto it = trial->coords_by_conformer_id.find(kv.first);
      ASSERT_NE(it, trial->coords_by_conformer_id.end()) << "Missing frame id=" << kv.first << " for t=" << t;
      const auto &a = kv.second;
      const auto &b = it->second;
      ASSERT_EQ(a.size(), b.size()) << "Vector size mismatch at frame id=" << kv.first << " for t=" << t;
      for (std::size_t i = 0; i < a.size(); ++i) {
        EXPECT_LE(sq_distance(a[i], b[i]), tol)
            << "Coordinate mismatch at frame id=" << kv.first << ", idx=" << i << " for t=" << t;
      }
    }
  }
}
