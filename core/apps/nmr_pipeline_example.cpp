#include <cassert>
#include <iomanip>
#include <iostream>
#include <logging.hpp>
#include <stdexcept>
#include <string>

#include <GraphMol/Conformer.h>

#include "lahuta.hpp"
#include "pipeline/ingestion.hpp"
#include "pipeline/pipeline_item.hpp"
#include "sources/nmr.hpp"
#include "topology.hpp"

int main() {
  using namespace lahuta;

  lahuta::Logger::get_instance().set_log_level(lahuta::Logger::LogLevel::Debug);
  const std::string session_id = "nmr/2LNL";
  IngestDescriptor descriptor{session_id, NMRRef{"core/data/2LNL.cif.gz"}};

  NMRRealizer realizer;

  std::shared_ptr<const Luni> cached_system;
  std::shared_ptr<const Topology> cached_topology;
  std::size_t frame_counter = 0;

  realizer.realize(descriptor, [&](PipelineItem item) {
    auto session = item.session;
    if (!session) {
      throw std::runtime_error("NMR realization did not attach a session");
    }

    // Each frame acquires a permit to respect per-session backpressure limits.
    [[maybe_unused]] auto permit = session->acquire_permit();

    if (!cached_system) {
      cached_system = session->get_or_load_system();
      std::cout << "System constructed once for session " << session->get_session_id() << " ("
                << cached_system.get() << ")\n";

      TopologyBuildingOptions opts{};
      cached_topology = session->get_or_load_topology(opts);
      std::cout << "Topology built once at address " << cached_topology.get() << "\n";
    } else {
      auto system_again = session->get_or_load_system();
      assert(system_again.get() == cached_system.get());
      auto topo_again = session->get_or_load_topology(TopologyBuildingOptions{});
      assert(topo_again.get() == cached_topology.get());
    }

    auto coordinates = item.frame->load_coordinates();
    std::shared_ptr<RDGeom::POINT3D_VECT> slab =
        coordinates.shared_positions
            ? std::const_pointer_cast<RDGeom::POINT3D_VECT>(coordinates.shared_positions)
            : std::make_shared<RDGeom::POINT3D_VECT>(std::move(coordinates.positions));
    RDKit::Conformer conformer;
    conformer.set3D(true);
    conformer.bindExternalPositions(std::move(slab));
    const RDKit::Conformer &cref = conformer;
    const auto &positions = cref.getPositions();
    assert(!positions.empty());

    const auto &first = positions.front();
    std::cout << "Frame " << std::setw(3) << item.conformer_id << ": first atom coordinates -> ("
              << std::fixed << std::setprecision(3) << first.x << ", " << first.y << ", " << first.z << ")"
              << std::endl;

    ++frame_counter;
  });

  std::cout << "Processed " << frame_counter << " frames from 2LNL" << std::endl;
  return 0;
}
