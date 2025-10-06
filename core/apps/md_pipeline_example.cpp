#include <cassert>
#include <iomanip>
#include <iostream>
#include <logging.hpp>
#include <stdexcept>
#include <string>

#include "pipeline/ingestion.hpp"
#include "pipeline/pipeline_item.hpp"
#include "sources/trajectory.hpp"
#include "topology.hpp"
#include <GraphMol/Conformer.h>

namespace {
constexpr const char *kStructurePath = "/Users/bsejdiu/Downloads/Atomistic_md_r1_cdl.tpr.gro";
constexpr const char *kXtcPath = "/Users/bsejdiu/Downloads/Atomistic_md_r1_cdl.pbc.xtc";
} // namespace

int main() {
  using namespace lahuta;

  lahuta::Logger::get_instance().set_log_level(lahuta::Logger::LogLevel::Debug);
  const std::string session_id = "md/example";
  IngestDescriptor descriptor{
      session_id,
      MDRef{kStructurePath, {kXtcPath}},
  };

  TrajectoryRealizer realizer;

  std::shared_ptr<const Luni> cached_system;
  std::shared_ptr<const Topology> cached_topology;
  std::size_t frame_counter = 0;

  while (frame_counter < 30) {
    auto start = std::chrono::high_resolution_clock::now();

    auto item = realizer.next(descriptor);
    if (!item) break;

    auto session = item->session;
    if (!session) {
      throw std::runtime_error("Trajectory realization did not attach a session");
    }

    [[maybe_unused]] auto permit = session->acquire_permit();

    if (!cached_system) {
      cached_system = session->get_or_load_system();
      std::cout << "System constructed once for session " << session->get_session_id() << " ("
                << cached_system.get() << ")\n";

      // TopologyBuildingOptions opts{};
      // cached_topology = session->get_or_load_topology(opts);
      // if (cached_topology) {
      //   std::cout << "Topology built once at address " << cached_topology.get() << "\n";
      // } else {
      //   std::cout << "Topology unavailable for structure-derived system\n";
      // }
    } else {
      // auto system_again = session->get_or_load_system();
      // assert(system_again.get() == cached_system.get());
      // auto topo_again = session->get_or_load_topology(TopologyBuildingOptions{});
      // if (cached_topology) {
      //   assert(topo_again.get() == cached_topology.get());
      // }
    }

    auto coordinates = item->frame->load_coordinates();
    std::shared_ptr<RDGeom::POINT3D_VECT> slab =
        coordinates.shared_positions
            ? std::const_pointer_cast<RDGeom::POINT3D_VECT>(coordinates.shared_positions)
            : std::make_shared<RDGeom::POINT3D_VECT>(std::move(coordinates.positions));
    RDKit::Conformer conformer;
    conformer.set3D(true);
    conformer.bindExternalPositions(std::move(slab));
    const RDKit::Conformer &cref = conformer; // const-view for read-only access
    const auto &positions = cref.getPositions();
    if (positions.empty()) {
      throw std::runtime_error("Trajectory frame produced no coordinates");
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> diff = end - start;
    std::cout << "Frame " << std::setw(6) << item->conformer_id << " loaded in " << std::fixed
              << std::setprecision(2) << diff.count() << " ms ";

    const auto &first = positions.front();
    std::cout << "Frame " << std::setw(6) << item->conformer_id << ": first atom => (" << std::fixed
              << std::setprecision(3) << first.x << ", " << first.y << ", " << first.z << ")";
    if (item->timestamp_ps) {
      std::cout << " time_ps=" << *item->timestamp_ps;
    }
    std::cout << std::endl;

    ++frame_counter;
  }

  std::cout << "Processed " << frame_counter << " frames from trajectory" << std::endl;
  return 0;
}
