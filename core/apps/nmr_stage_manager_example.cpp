#include <analysis/contacts/computation.hpp>
#include <io/sinks/ndjson.hpp>
#include <iostream>
#include <memory>
#include <mutex>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

#include <GraphMol/Conformer.h>

#include "analysis/topology/computation.hpp"
#include "compute/dependency.hpp"
#include "compute/result.hpp"
#include "lahuta.hpp"
#include "pipeline/compute/context.hpp"
#include "pipeline/dynamic/manager.hpp"
#include "pipeline/ingestion.hpp"
#include "pipeline/stream_session.hpp"
#include "topology.hpp"

// clang-format off
namespace {

using namespace lahuta;
using lahuta::IngestDescriptor;
using namespace lahuta::topology::compute;
using namespace lahuta::pipeline::compute;
using lahuta::sources::IDescriptor;

class SingleNMRDescriptor final : public IDescriptor {
public:
  explicit SingleNMRDescriptor(std::string path) : path_(std::move(path)) {}

  std::optional<IngestDescriptor> next() override {
    if (done_) return std::nullopt;
    done_ = true;
    IngestDescriptor desc;
    desc.id = path_;
    desc.origin = NMRRef{path_};
    return desc;
  }

  void reset() override { done_ = false; }

private:
  std::string path_;
  bool done_ = false;
};

struct FrameSummaryParams : ParameterBase<FrameSummaryParams> {
  static constexpr ParameterInterface::TypeId TYPE_ID = 200;
};

class FrameSummaryComputation final : public ReadWriteComputation<PipelineContext, FrameSummaryParams, FrameSummaryComputation> {
public:
  using Base = ReadWriteComputation<PipelineContext, FrameSummaryParams, FrameSummaryComputation>;
  FrameSummaryComputation() : Base(FrameSummaryParams{}) {}

  static constexpr ComputationLabel label{"frame_summary"};
  using dependencies = Dependencies<Dependency<analysis::topology::BuildTopologyComputation, void>>;

  ComputationResult execute_typed(DataContext<PipelineContext, Mut::ReadWrite> &context, const FrameSummaryParams &) {
    auto &data = context.data();
    auto *task_ctx = data.ctx;

    auto system = task_ctx ? task_ctx->get_object<const Luni>("system") : nullptr;
    if (!system && data.session) {
      system = data.session->get_or_load_system();
      if (task_ctx && system) {
        task_ctx->set_object<const Luni>("system", system);
      }
    }
    if (!system) {
      return ComputationResult(ComputationError("FrameSummaryComputation requires a system"));
    }

    auto topology_ptr = task_ctx ? task_ctx->get_object<const Topology>("topology") : nullptr;
    if (!topology_ptr && data.session) {
      TopologyBuildingOptions opts{};
      topology_ptr = data.session->get_or_load_topology(opts);
      if (task_ctx && topology_ptr) {
        task_ctx->set_object<const Topology>("topology", topology_ptr);
      }
    }

    bool same_system = !first_system_ || first_system_.get() == system.get();
    if (!first_system_) first_system_ = system;
    bool same_topology = true;
    if (topology_ptr) {
      if (!first_topology_) {
        first_topology_ = topology_ptr;
      } else {
        same_topology = first_topology_.get() == topology_ptr.get();
      }
    } else {
      same_topology = false;
    }

    RDGeom::Point3D first_atom{};
    if (data.frame && data.session) {
      auto coords = data.frame->load_coordinates();
      std::shared_ptr<RDGeom::POINT3D_VECT> slab =
          coords.shared_positions ? std::const_pointer_cast<RDGeom::POINT3D_VECT>(coords.shared_positions)
                                  : std::make_shared<RDGeom::POINT3D_VECT>(std::move(coords.positions));
      RDKit::Conformer conf;
      conf.set3D(true);
      conf.bindExternalPositions(std::move(slab));
      const RDKit::Conformer &cref = conf;
      first_atom = cref.getAtomPos(0);
    } else {
      const auto &conf = system->get_conformer();
      first_atom = conf.getAtomPos(0);
    }

    auto frame_meta = task_ctx ? task_ctx->get_object<FrameMetadata>("lahuta.frame") : nullptr;

    static std::once_flag header_once;
    std::call_once(header_once, []() {
      std::cout << "StageManager-driven NMR streaming example" << std::endl;
      std::cout << "System/topology are built once per session. Coordinates update per frame." << std::endl;
    });

    std::ostringstream oss;
    oss << "Frame " << data.conformer_id;
    if (frame_meta) {
      oss << " session='" << frame_meta->session_id << "'";
      if (frame_meta->timestamp_ps) {
        oss << " time_ps=" << *frame_meta->timestamp_ps;
      }
    } else {
      oss << " (single-structure)";
    }
    oss << " system_ptr=" << static_cast<const void *>(system.get());
    oss << (same_system ? " (reused)" : " (!! new instance)");
    if (topology_ptr) {
      oss << " topology_ptr=" << static_cast<const void *>(topology_ptr.get());
      oss << (same_topology ? " (reused)" : " (!! new instance)");
    } else {
      oss << " topology_ptr=<none>";
    }
    oss << " first_atom=(" << first_atom.x << ", " << first_atom.y << ", " << first_atom.z << ")";

    std::cout << oss.str() << std::endl;
    return ComputationResult(true);
  }

private:
  std::shared_ptr<const Luni> first_system_;
  std::shared_ptr<const Topology> first_topology_;
};

} // namespace

int main(int argc, char **argv) {
  Logger::get_instance().set_log_level(Logger::LogLevel::Debug);

  std::string input = argc > 1 ? argv[1] : "core/data/2LNL.cif.gz";

  auto source = std::make_unique<SingleNMRDescriptor>(input);
  pipeline::dynamic::StageManager manager(std::move(source));
  manager.set_auto_builtins(true);

  ContactsParams p{};
  manager.add_computation(
      "contacts",
      {},
      [label = std::string("contacts"), p]() {
        return std::make_unique<analysis::contacts::ContactsComputation>(label, p);
      },
      /*thread_safe=*/true);

  auto contacts_data_sink = std::make_shared<lahuta::pipeline::dynamic::NdjsonFileSink>("contacts_data.json");
  manager.connect_sink("contacts", contacts_data_sink);

  manager.add_computation("frame_summary", {"topology"}, [] {
    return std::make_unique<FrameSummaryComputation>();
  });

  manager.run(16);
  return 0;
}
