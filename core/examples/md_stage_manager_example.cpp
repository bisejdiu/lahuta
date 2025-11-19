#include <analysis/contacts/computation.hpp>
#include <analysis/contacts/hooks.hpp>
#include <iostream>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <vector>

#include <GraphMol/Conformer.h>

#include "analysis/system/computation.hpp"
#include "compute/dependency.hpp"
#include "compute/result.hpp"
#include "pipeline/compute/context.hpp"
#include "pipeline/dynamic/manager.hpp"
#include "pipeline/ingestion.hpp"
#include "pipeline/stream_session.hpp"
#include "sinks/ndjson.hpp"

// clang-format off
namespace {

using namespace lahuta;
using lahuta::sources::IDescriptor;
using namespace lahuta::pipeline::dynamic;
using namespace lahuta::pipeline::compute;

/*constexpr const char* kStructurePath = "/Users/bsejdiu/Downloads/Atomistic_md_r1_cdl.tpr.gro";*/
/*constexpr const char* kXtcPath = "/Users/bsejdiu/Downloads/Atomistic_md_r1_cdl.pbc.xtc";*/
constexpr const char *kStructurePath = "/Users/bsejdiu/projects/lahuta_dev/lahuta/gpcrmd/10828_dyn_85.pdb";
constexpr const char *kXtcPath = "/Users/bsejdiu/projects/lahuta_dev/lahuta/gpcrmd/10824_trj_85.xtc";

class SingleTrajectoryDescriptor final : public IDescriptor {
public:
  SingleTrajectoryDescriptor(std::string structure, std::vector<std::string> xtcs)
      : structure_(std::move(structure)), xtcs_(std::move(xtcs)) {}

  std::optional<IngestDescriptor> next() override {
    if (done_) return std::nullopt;
    done_ = true;
    IngestDescriptor desc;
    desc.id = structure_;
    desc.origin = MDRef{structure_, xtcs_};
    return desc;
  }

  void reset() override { done_ = false; }

private:
  std::string structure_;
  std::vector<std::string> xtcs_;
  bool done_ = false;
};

struct TrajectorySummaryParams : ParameterBase<TrajectorySummaryParams> {
  static constexpr ParameterInterface::TypeId TYPE_ID = 201;
  static constexpr const char* DEFAULT_CHANNEL = "trajectory_summary";
  std::string channel = DEFAULT_CHANNEL;
};

class TrajectorySummaryComputation final : public ReadWriteComputation<PipelineContext, TrajectorySummaryParams, TrajectorySummaryComputation> {
public:
  using Base = ReadWriteComputation<PipelineContext, TrajectorySummaryParams, TrajectorySummaryComputation>;
  TrajectorySummaryComputation() : Base(TrajectorySummaryParams{}) {}

  static constexpr ComputationLabel label{TrajectorySummaryParams::DEFAULT_CHANNEL};
  using dependencies = Dependencies<Dependency<analysis::system::SystemReadComputation, void>>;

  ComputationResult execute_typed(DataContext<PipelineContext, Mut::ReadWrite> &context, const TrajectorySummaryParams &p) {
    auto &data = context.data();
    auto *task_ctx = data.ctx;

    std::shared_ptr<const Luni> system = task_ctx ? task_ctx->system() : nullptr;
    if (!system) {
      return ComputationResult(ComputationError("TrajectorySummaryComputation requires a system"));
    }

    std::vector<RDGeom::Point3D> first_5_atoms{};
    if (data.frame && data.session) {
      auto coords = data.frame->load_coordinates();
      std::shared_ptr<RDGeom::POINT3D_VECT> slab =
          coords.shared_positions ? std::const_pointer_cast<RDGeom::POINT3D_VECT>(coords.shared_positions)
                                  : std::make_shared<RDGeom::POINT3D_VECT>(std::move(coords.positions));
      RDKit::Conformer conf;
      conf.set3D(true);
      conf.bindExternalPositions(std::move(slab));
      const RDKit::Conformer &cref = conf; // const-view to avoid mutating accessors
      const int num_atoms = std::min(5, static_cast<int>(cref.getNumAtoms()));
      for (int i = 0; i < num_atoms; ++i) {
        first_5_atoms.push_back(cref.getAtomPos(i));
      }
    } else {
      const auto &conf = system->get_conformer();
      const int num_atoms = std::min(5, static_cast<int>(conf.getNumAtoms()));
      for (int i = 0; i < num_atoms; ++i) {
        first_5_atoms.push_back(conf.getAtomPos(i));
      }
    }

    auto frame_meta = task_ctx ? task_ctx->frame_metadata() : nullptr;

    static std::once_flag header_once;
    std::call_once(header_once, []() {
      std::cout << "StageManager driven MD streaming example" << std::endl;
      std::cout << "System/topology are reused per session. Coordinates update for each frame." << std::endl;
    });

    pipeline::dynamic::EmissionList out;
    auto json = lahuta::analysis::contacts::build_contacts_summary_json(*data.ctx);
    out.push_back({p.channel, std::move(json)});

    return ComputationResult(std::move(out));
  }
};

} // namespace

int main() {
  Logger::get_instance().set_log_level(Logger::LogLevel::Warn);
  auto source = std::make_unique<SingleTrajectoryDescriptor>(kStructurePath, std::vector<std::string>{kXtcPath});
  StageManager manager(std::move(source));
  manager.set_auto_builtins(true);

  // manager.add_computation(
  //     "system",
  //     {}, // No dependencies
  //     []() {
  //       SystemReadParams sys_params_;
  //       return std::make_unique<analysis::system::SystemReadComputation>(sys_params_);
  //     },
  //     /*thread_safe=*/true
  // );

  // manager.add_computation(
  //     "topology",
  //     {"system"}, // No dependencies
  //     []() {
  //       BuildTopologyParams top_params;
  //       return std::make_unique<analysis::topology::BuildTopologyComputation>(top_params);
  //     },
  //     /*thread_safe=*/true
  // );

  ContactsParams p{};
  manager.add_computation(
      "contacts",
      {},
      [label = std::string("contacts"), p]() {
        return std::make_unique<analysis::contacts::ContactsComputation>(label, p);
      },
      /*thread_safe=*/true);

  manager.add_computation("trajectory_summary", {"contacts"}, [] {
    return std::make_unique<TrajectorySummaryComputation>();
  });

  auto contacts_data_sink = std::make_shared<lahuta::pipeline::dynamic::NdjsonFileSink>("contacts_data.json");
  manager.connect_sink("contacts", contacts_data_sink);

  auto contacts_summary_sink = std::make_shared<lahuta::pipeline::dynamic::NdjsonFileSink>("contacts_summary.json");
  manager.connect_sink(TrajectorySummaryParams::DEFAULT_CHANNEL, contacts_summary_sink);

  manager.run(16);
  return 0;
}
