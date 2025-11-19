#ifndef LAHUTA_COMPUTE_TOPOLOGY_SNAPSHOT_HPP
#define LAHUTA_COMPUTE_TOPOLOGY_SNAPSHOT_HPP

#include <optional>
#include <stdexcept>

#include <GraphMol/Conformer.h>

#include "pipeline/dynamic/types.hpp"
#include "topology.hpp"

// clang-format off
namespace lahuta::compute {

// Immutable view pairing a Topology with the active coordinates
struct TopologySnapshot {
  const Topology& topo;
  const RDKit::Conformer& conf;
};

// Construct from a conformer
inline TopologySnapshot snapshot_of(const Topology& topo, const RDKit::Conformer& conf) {
  return TopologySnapshot{topo, conf};
}

// Construct from a Topology and optional TaskContext
inline TopologySnapshot snapshot_of(const Topology& topo, const pipeline::dynamic::TaskContext* tctx = nullptr) {
  if (tctx) {
    if (auto conf = tctx->conformer()) {
      return TopologySnapshot{topo, *conf};
    }
  }
  return TopologySnapshot{topo, topo.conformer()};
}

// Non-null TaskContext reference
inline TopologySnapshot snapshot_of(const Topology& topo, const pipeline::dynamic::TaskContext& tctx) {
  if (auto conf = tctx.conformer()) {
    return TopologySnapshot{topo, *conf};
  }
  return TopologySnapshot{topo, topo.conformer()};
}

// Return a snapshot if Topology exists in context, otherwise std::nullopt.
inline std::optional<TopologySnapshot> try_topology_snapshot(const pipeline::dynamic::TaskContext& tctx) {
  auto top = tctx.topology();
  if (!top) return std::nullopt;
  if (auto conf = tctx.conformer()) {
    return TopologySnapshot{*top, *conf};
  }
  return TopologySnapshot{*top, top->conformer()};
}

// Require a snapshot from context, throws if Topology is missing.
inline TopologySnapshot require_topology_snapshot(const pipeline::dynamic::TaskContext& tctx) {
  auto top = tctx.topology();
  if (!top) throw std::runtime_error("Topology not found in TaskContext");
  if (auto conf = tctx.conformer()) {
    return TopologySnapshot{*top, *conf};
  }
  return TopologySnapshot{*top, top->conformer()};
}

} // namespace lahuta::compute

#endif // LAHUTA_COMPUTE_TOPOLOGY_SNAPSHOT_HPP
