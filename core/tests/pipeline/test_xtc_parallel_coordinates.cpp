#include <atomic>
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
#include "logging.hpp"
#include "pipeline/compute/context.hpp"
#include "pipeline/dynamic/manager.hpp"
#include "pipeline/ingestion.hpp"
#include "sources/descriptor.hpp"

using namespace lahuta;
using namespace lahuta::topology::compute;
using namespace lahuta::pipeline::compute;

// clang-format off
namespace {

constexpr const char *gpcrmd_dir = "/Users/bsejdiu/projects/lahuta_dev/lahuta/gpcrmd";
std::string pdb_path = std::string(gpcrmd_dir) + "/10828_dyn_85.pdb";
std::string xtc_path = std::string(gpcrmd_dir) + "/10824_trj_85.xtc";

struct SingleTrajectoryDescriptor final : sources::IDescriptor {
  std::string structure;
  std::vector<std::string> xtcs;
  bool done{false};

  SingleTrajectoryDescriptor(std::string s, std::vector<std::string> x)
      : structure(std::move(s)), xtcs(std::move(x)) {}

  std::optional<IngestDescriptor> next() override {
    if (done) return std::nullopt;
    done = true;
    IngestDescriptor d;
    d.id = structure;
    d.origin = MDRef{structure, xtcs};
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
  bool indices_ready{false};
  std::vector<unsigned> selected_indices; // size <= 100
  std::unordered_map<std::size_t, std::vector<RDGeom::Point3D>> coords_by_conformer_id;
  std::atomic<size_t> frames{0};

  void reset() {
    std::scoped_lock lk(mtx);
    indices_ready = false;
    selected_indices.clear();
    coords_by_conformer_id.clear();
    frames = 0;
  }
};

static std::vector<unsigned> select_random_unique_indices(std::size_t num_atoms, std::size_t k, uint64_t seed) {
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

struct CoordsCaptureParams : ParameterBase<CoordsCaptureParams> {
  static constexpr ParameterInterface::TypeId TYPE_ID = 204; // local test id
};

class CoordsCaptureComputation final : public ReadWriteComputation<PipelineContext, CoordsCaptureParams, CoordsCaptureComputation> {
public:
  using Base = ReadWriteComputation<PipelineContext, CoordsCaptureParams, CoordsCaptureComputation>;
  explicit CoordsCaptureComputation(std::shared_ptr<CoordinateCollector> collector)
      : Base(CoordsCaptureParams{}), collector_(std::move(collector)) {}

  static constexpr ComputationLabel label{"coords_capture"};
  using dependencies = Dependencies<Dependency<analysis::system::SystemReadComputation, void>>;

  ComputationResult execute_typed(DataContext<PipelineContext, Mut::ReadWrite> &ctx, const CoordsCaptureParams &) {
    auto &data = ctx.data();

    if (!data.frame) return ComputationResult(true);

    // Build a conformer view bound to this frame's coordinates
    auto view = data.frame->load_coordinates();
    auto slab = view.shared_positions ? std::const_pointer_cast<RDGeom::POINT3D_VECT>(view.shared_positions)
                                      : std::make_shared<RDGeom::POINT3D_VECT>(std::move(view.positions));
    RDKit::Conformer conf;
    conf.set3D(true);
    conf.bindExternalPositions(std::move(slab));
    const RDKit::Conformer &cref = conf; // const-view to avoid mutating accessors

    if (!collector_->indices_ready) {
      std::scoped_lock lk(collector_->mtx);
      if (!collector_->indices_ready) {
        const std::size_t natoms = static_cast<std::size_t>(cref.getNumAtoms());
        collector_->selected_indices = select_random_unique_indices(natoms, 100, /*seed=*/1337);
        collector_->indices_ready = true;
      }
    }

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
    return ComputationResult(true);
  }

private:
  std::shared_ptr<CoordinateCollector> collector_;
};

} // namespace

TEST(ParallelXtcCoordinates, ConsistencyAcrossThreads) {
  Logger::get_instance().set_log_level(Logger::LogLevel::Info);

  // Oracle
  auto baseline = std::make_shared<CoordinateCollector>();
  baseline->reset();
  {
    auto src = std::make_unique<SingleTrajectoryDescriptor>(std::string(pdb_path), std::vector<std::string>{std::string(xtc_path)});
    pipeline::dynamic::StageManager mgr(std::move(src));
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

    auto src = std::make_unique<SingleTrajectoryDescriptor>(std::string(pdb_path), std::vector<std::string>{std::string(xtc_path)});
    pipeline::dynamic::StageManager mgr(std::move(src));
    mgr.set_auto_builtins(true);
    mgr.add_computation("coords_capture", {"system"}, [trial] {
      return std::make_unique<CoordsCaptureComputation>(trial);
    });
    mgr.run(static_cast<std::size_t>(t));

    // Basic invariants
    ASSERT_EQ(trial->frames.load(), baseline_frames) << "Thread count " << t << " produced different frame count";
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
        EXPECT_LE(sq_distance(a[i], b[i]), tol) << "Coordinate mismatch at frame id=" << kv.first
                                                << ", idx=" << i << " for t=" << t;
      }
    }
  }
}
