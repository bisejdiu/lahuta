#include <atomic>
#include <cmath>
#include <filesystem>
#include <mutex>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <GraphMol/Conformer.h>
#include <gtest/gtest.h>

#include "analysis/topology/computation.hpp"
#include "compute/dependency.hpp"
#include "compute/result.hpp"
#include "entities/records.hpp"
#include "logging.hpp"
#include "pipeline/compute/context.hpp"
#include "pipeline/dynamic/keys.hpp"
#include "pipeline/dynamic/manager.hpp"
#include "pipeline/ingestion.hpp"
#include "sources/descriptor.hpp"

using namespace lahuta;
using namespace lahuta::topology::compute;
using namespace lahuta::pipeline::compute;

namespace {

namespace fs = std::filesystem;

static fs::path locate_simulation_file(const std::string &filename) {
#ifdef LAHUTA_SIMDB_DIR
  {
    fs::path from_cmake{LAHUTA_SIMDB_DIR};
    fs::path cand = from_cmake / filename;
    if (fs::exists(cand)) return cand;
  }
#endif
  fs::path p = fs::current_path();
  for (int i = 0; i < 6; ++i) {
    fs::path cand_core = p / "core" / "data" / "simulationdatabase" / filename;
    if (fs::exists(cand_core)) return cand_core;

    if (!p.has_parent_path()) break;
    p = p.parent_path();
  }
  return {};
}

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

// Oracle capturing stable topology containers
struct TopologyOracle {
  struct ResidueRec {
    std::string chain;
    int number;
    std::string name;
    std::string alt;
    std::vector<unsigned> atom_idxs;
  };
  struct RingRecO {
    std::vector<unsigned> atom_idxs;
    bool aromatic;
  };
  struct GroupRecO {
    uint32_t a_type;
    int type;
    std::vector<unsigned> atom_idxs;
  };

  std::mutex mtx;
  bool initialized{false};
  std::vector<ResidueRec> residues;
  std::vector<RingRecO> rings;
  std::vector<GroupRecO> groups;
  std::atomic<size_t> frames{0};
  std::atomic<size_t> errors{0};
  std::atomic<size_t> changed_frames{0};
  // baseline ring centers for motion detection across frames
  bool centers_ready{false};
  std::vector<RDGeom::Point3D> last_ring_centers;
  bool positions_ready{false};
  std::vector<RDGeom::Point3D> last_positions;
  bool frame_index_ready{false};
  std::size_t last_frame_index{0};
  // topology object identity validation
  std::atomic<const Topology *> last_topology_ptr{nullptr};

  void reset() {
    std::scoped_lock lk(mtx);
    initialized = false;
    centers_ready = false;
    positions_ready = false;
    frame_index_ready = false;
    residues.clear();
    rings.clear();
    groups.clear();
    last_ring_centers.clear();
    last_positions.clear();
    last_frame_index = 0;
    frames.store(0);
    errors.store(0);
    changed_frames.store(0);
    last_topology_ptr.store(nullptr);
  }
} g_oracle;

static std::vector<unsigned>
atoms_to_indices(const std::vector<std::reference_wrapper<const RDKit::Atom>> &atoms) {
  std::vector<unsigned> v;
  v.reserve(atoms.size());
  for (const auto &a : atoms)
    v.push_back(static_cast<unsigned>(a.get().getIdx()));
  return v;
}

static std::vector<unsigned> residue_atoms_to_indices(const std::vector<const RDKit::Atom *> &atoms) {
  std::vector<unsigned> v;
  v.reserve(atoms.size());
  for (const auto *a : atoms)
    v.push_back(static_cast<unsigned>(a->getIdx()));
  return v;
}

// Direct geometric recompute from a conformer
static RDGeom::Point3D direct_center(const RDKit::Conformer &conf, const std::vector<unsigned> &idxs) {
  RDGeom::Point3D c{0.0, 0.0, 0.0};
  if (idxs.empty()) return c;
  for (auto idx : idxs)
    c += conf.getAtomPos(idx);
  c /= static_cast<double>(idxs.size());
  return c;
}

static RDGeom::Point3D direct_normal_first3(const RDKit::Conformer &conf, const std::vector<unsigned> &idxs) {
  RDGeom::Point3D n{0.0, 0.0, 0.0};
  if (idxs.size() < 3) return n;
  auto p0 = conf.getAtomPos(idxs[0]);
  auto p1 = conf.getAtomPos(idxs[1]);
  auto p2 = conf.getAtomPos(idxs[2]);
  n = (p1 - p0).crossProduct(p2 - p0);
  const double len = std::sqrt(n.lengthSq());
  if (len > 0.0) n /= len;
  return n;
}

static double sq_distance(const RDGeom::Point3D &a, const RDGeom::Point3D &b) {
  double dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
  return dx * dx + dy * dy + dz * dz;
}

struct OracleCheckParams : ParameterBase<OracleCheckParams> {
  static constexpr ParameterInterface::TypeId TYPE_ID = 203; // local test id (fits in uint8)
};

class OracleCheckComputation final
    : public ReadWriteComputation<PipelineContext, OracleCheckParams, OracleCheckComputation> {
public:
  using Base = ReadWriteComputation<PipelineContext, OracleCheckParams, OracleCheckComputation>;
  OracleCheckComputation() : Base(OracleCheckParams{}) {}

  static constexpr ComputationLabel label{"oracle_check"};
  using dependencies = Dependencies<Dependency<analysis::topology::BuildTopologyComputation, void>>;

  ComputationResult
  execute_typed(DataContext<PipelineContext, Mut::ReadWrite> &ctx, const OracleCheckParams &) {
    auto &data = ctx.data();
    try {
      // system and topology
      auto task_ctx = data.ctx;
      auto topo = task_ctx ? task_ctx->topology() : nullptr;
      if (!topo && data.session) {
        TopologyBuildingOptions opts{};
        topo = data.session->get_or_load_topology(opts);
      }
      if (!topo) return ComputationResult(ComputationError("oracle_check requires topology"));

      // Validate topology object identity
      {
        auto current_topo_ptr = topo.get();
        auto expected_ptr = g_oracle.last_topology_ptr.load();

        if (expected_ptr == nullptr) {
          // First frame - store the topology pointer
          g_oracle.last_topology_ptr.store(current_topo_ptr);
        } else if (current_topo_ptr != expected_ptr) {
          // Topology object changed across frames, caching failed
          ++g_oracle.errors;
          return ComputationResult(ComputationError("Topology object identity changed across frames"));
        }
      }

      // conformer view
      auto conf_ptr = task_ctx ? task_ctx->conformer() : nullptr;
      RDKit::Conformer conf_tmp; // fallback
      if (!conf_ptr && data.frame) {
        auto view = data.frame->load_coordinates();
        std::shared_ptr<RDGeom::POINT3D_VECT> slab =
            view.shared_positions ? std::const_pointer_cast<RDGeom::POINT3D_VECT>(view.shared_positions)
                                  : std::make_shared<RDGeom::POINT3D_VECT>(std::move(view.positions));
        conf_tmp.set3D(true);
        conf_tmp.bindExternalPositions(std::move(slab));
      }
      const RDKit::Conformer &conf = conf_ptr ? *conf_ptr : conf_tmp;

      bool frame_moved = false;
      //
      // I track a few (maybe redundant) signals:
      //   - ring centers
      //   - overall atom positions
      //   - the frame indices themselves
      // Only when any of these change do I consider the frame to have "moved".
      //

      // Initialize oracle
      {
        std::scoped_lock lk(g_oracle.mtx);
        if (!g_oracle.initialized) {
          // residues
          const auto &res = topo->get_residues().get_residues();
          g_oracle.residues.reserve(res.size());
          for (const auto &r : res) {
            TopologyOracle::ResidueRec
                rr{r.chain_id, r.number, r.name, r.alt_loc, residue_atoms_to_indices(r.atoms)};
            g_oracle.residues.push_back(std::move(rr));
          }
          // rings
          const auto &ring_recs = topo->records<RingRec>();
          g_oracle.rings.reserve(ring_recs.size());
          for (const auto &r : ring_recs) {
            g_oracle.rings.push_back(TopologyOracle::RingRecO{atoms_to_indices(r.atoms), r.aromatic});
          }
          // groups
          const auto &group_recs = topo->records<GroupRec>();
          g_oracle.groups.reserve(group_recs.size());
          for (const auto &g : group_recs) {
            g_oracle.groups.push_back(
                TopologyOracle::GroupRecO{
                    static_cast<uint32_t>(g.a_type),
                    static_cast<int>(g.type),
                    atoms_to_indices(g.atoms)});
          }
          g_oracle.initialized = true;
        }
      }

      // verify stability: residues
      {
        const auto &res = topo->get_residues().get_residues();
        if (res.size() != g_oracle.residues.size())
          ++g_oracle.errors;
        else {
          for (size_t i = 0; i < res.size(); ++i) {
            const auto &r = res[i];
            const auto &o = g_oracle.residues[i];
            if (r.chain_id != o.chain || r.number != o.number || r.name != o.name || r.alt_loc != o.alt) {
              ++g_oracle.errors;
              break;
            }
            if (residue_atoms_to_indices(r.atoms) != o.atom_idxs) {
              ++g_oracle.errors;
              break;
            }
          }
        }
      }

      // verify stability: rings
      {
        const auto &ring_recs = topo->records<RingRec>();
        if (ring_recs.size() != g_oracle.rings.size())
          ++g_oracle.errors;
        else {
          double method_delta = 0.0; // difference between accessor vs direct recompute (same-frame)
          std::vector<RDGeom::Point3D> current_centers;
          current_centers.reserve(ring_recs.size());
          for (size_t i = 0; i < ring_recs.size(); ++i) {
            const auto &r = ring_recs[i];
            const auto &o = g_oracle.rings[i];
            auto idxs = atoms_to_indices(r.atoms);
            if (idxs != o.atom_idxs || r.aromatic != o.aromatic) {
              ++g_oracle.errors;
              break;
            }
            // recompute geometry
            auto c1 = r.center(conf);
            auto n1 = r.normal(conf);
            auto c2 = direct_center(conf, idxs);
            auto n2 = direct_normal_first3(conf, idxs);
            const double tol = 1e-6;
            if (sq_distance(c1, c2) > 1e-4) {
              ++g_oracle.errors;
              break;
            }
            double dp = n1.dotProduct(n2);
            if (std::abs(dp) < 1.0 - 1e-4) {
              ++g_oracle.errors;
              break;
            }
            method_delta += sq_distance(c1, c2);
            current_centers.push_back(c1);
          }
          // motion detection across frames: compare to previous centers
          {
            std::scoped_lock lk(g_oracle.mtx);
            if (!g_oracle.centers_ready) {
              g_oracle.last_ring_centers = std::move(current_centers);
              g_oracle.centers_ready = true;
            } else if (current_centers.size() == g_oracle.last_ring_centers.size()) {
              const double move_tol = 1e-12;
              bool ring_moved = false;
              for (size_t i = 0; i < current_centers.size(); ++i) {
                if (sq_distance(current_centers[i], g_oracle.last_ring_centers[i]) > move_tol) {
                  ring_moved = true;
                  break;
                }
              }
              if (ring_moved) frame_moved = true;
              g_oracle.last_ring_centers.swap(current_centers);
            }
          }
        }
      }

      // verify stability: groups
      {
        const auto &group_recs = topo->records<GroupRec>();
        if (group_recs.size() != g_oracle.groups.size())
          ++g_oracle.errors;
        else {
          for (size_t i = 0; i < group_recs.size(); ++i) {
            const auto &g = group_recs[i];
            const auto &o = g_oracle.groups[i];
            auto idxs = atoms_to_indices(g.atoms);
            if (static_cast<uint32_t>(g.a_type) != o.a_type || static_cast<int>(g.type) != o.type
                || idxs != o.atom_idxs) {
              ++g_oracle.errors;
              break;
            }
            auto c1 = g.center(conf);
            auto c2 = direct_center(conf, idxs);
            if (sq_distance(c1, c2) > 1e-4) {
              ++g_oracle.errors;
              break;
            }
          }
        }
      }

      if (data.frame) {
        const auto frame_idx = data.frame->index();
        std::scoped_lock lk(g_oracle.mtx);
        if (!g_oracle.frame_index_ready) {
          g_oracle.last_frame_index = frame_idx;
          g_oracle.frame_index_ready = true;
        } else if (g_oracle.last_frame_index != frame_idx) {
          frame_moved = true;
          g_oracle.last_frame_index = frame_idx;
        }
      }

      // capture overall coordinate motion as a fallback
      {
        const auto &positions = conf.getPositions();
        std::vector<RDGeom::Point3D> current_positions;
        current_positions.reserve(positions.size());
        for (const auto &p : positions)
          current_positions.push_back(p);

        std::scoped_lock lk(g_oracle.mtx);
        if (!g_oracle.positions_ready) {
          g_oracle.last_positions = std::move(current_positions);
          g_oracle.positions_ready = true;
        } else if (g_oracle.last_positions.size() == current_positions.size()) {
          const double pos_tol = 1e-12;
          bool coordinates_moved = false;
          for (size_t i = 0; i < current_positions.size(); ++i) {
            if (sq_distance(current_positions[i], g_oracle.last_positions[i]) > pos_tol) {
              coordinates_moved = true;
              break;
            }
          }
          if (coordinates_moved) frame_moved = true;
          g_oracle.last_positions.swap(current_positions);
        } else {
          ++g_oracle.errors;
        }
      }

      if (frame_moved) {
        g_oracle.changed_frames.fetch_add(1);
      }

      ++g_oracle.frames;
      return ComputationResult(true);
    } catch (const std::exception &e) {
      return ComputationResult(ComputationError(std::string("oracle_check exception: ") + e.what()));
    }
  }
};

} // namespace

TEST(TopologyCacheValidation, GPCRMD_Trajectory_Threads) {
  Logger::get_instance().set_log_level(Logger::LogLevel::Info);

  const auto structure_path  = locate_simulation_file("lysozyme.gro");
  const auto trajectory_path = locate_simulation_file("lysozyme.xtc");
  if (structure_path.empty() || trajectory_path.empty()) {
    GTEST_SKIP() << "core/data/simulationdatabase lysozyme data not available";
  }

  const std::vector<int> threads = {1, 4, 8, 16};
  for (int t : threads) {
    g_oracle.reset();

    auto src = std::make_unique<SingleTrajectoryDescriptor>(
        structure_path.string(),
        std::vector<std::string>{trajectory_path.string()});
    pipeline::dynamic::StageManager mgr(std::move(src));
    mgr.set_auto_builtins(true);
    mgr.add_computation("oracle_check", {"topology"}, [] {
      return std::make_unique<OracleCheckComputation>();
    });

    mgr.run(static_cast<std::size_t>(t));

    EXPECT_GT(g_oracle.frames.load(), 0u) << "No frames processed";
    EXPECT_EQ(g_oracle.errors.load(), 0u) << "Topology stability or geometry recompute failed";
    EXPECT_GT(g_oracle.changed_frames.load(), 0u) << "Expected motion in multi-frame trajectory";
  }
}
