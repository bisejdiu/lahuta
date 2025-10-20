#ifndef LAHUTA_ANALYSIS_TOPOLOGY_KERNEL_HPP
#define LAHUTA_ANALYSIS_TOPOLOGY_KERNEL_HPP

#include <memory>
#include <string>

#include "compute/result.hpp"
#include "lahuta.hpp"
#include "pipeline/compute/context.hpp"
#include "pipeline/compute/parameters.hpp"
#include "topology.hpp"

// clang-format off
namespace lahuta::analysis::topology {
using namespace pipeline::compute;

struct BuildTopologyKernel {
  static ComputationResult execute(DataContext<PipelineContext, Mut::ReadWrite>& context, const BuildTopologyParams& p) {
    try {
      auto& data = context.data();
      std::shared_ptr<const Luni> sys;
      if (data.ctx) sys = data.ctx->get_object<const Luni>(pipeline::CTX_SYSTEM_KEY);
      if (!sys && data.session) {
        auto shared = data.session->get_or_load_system();
        if (shared) {
          sys = shared;
          if (data.ctx) data.ctx->set_object<const Luni>(pipeline::CTX_SYSTEM_KEY, shared);
        }
      }
      if (!sys) return ComputationResult(ComputationError("BuildTopology requires system in context"));

      sys->set_atom_typing_method(p.atom_typing_method);

      bool is_model_mode = sys->is_model_origin();

      if (is_model_mode) {
        // Ensure topology built via model pathway if not already built
        TopologyBuildingOptions opts{};
        opts.mode = TopologyBuildMode::Model;
        opts.atom_typing_method = p.atom_typing_method;
        (void)sys->build_topology(opts); // idempotent if already built
        sys->enable_only(p.flags);
        auto topo = sys->get_topology();
        if (data.ctx && topo) data.ctx->set_object<const Topology>(pipeline::CTX_TOPOLOGY_KEY, topo);
        return ComputationResult(true);
      }

      // exec full topology comp
      sys->enable_only(p.flags);

      // edge case
      // allow disabling all computations while still materializing a Topology object in context.
      if (p.flags == TopologyComputation::None) {
        auto topo0 = sys->get_topology();
        if (!topo0) {
          sys->enable_only(TopologyComputation::None);
          topo0 = sys->get_topology();
        }
        if (data.ctx && topo0) data.ctx->set_object<const Topology>(pipeline::CTX_TOPOLOGY_KEY, topo0);
        return ComputationResult(true);
      }

      // normal full build
      TopologyBuildingOptions opts{};
      opts.compute_nonstandard_bonds = has_flag(p.flags, TopologyComputation::NonStandardBonds);
      opts.atom_typing_method = p.atom_typing_method;
      if (!sys->build_topology(opts)) return ComputationResult(ComputationError("BuildTopology failed"));
      sys->enable_only(p.flags);
      auto topo = sys->get_topology();
      if (data.ctx && topo) data.ctx->set_object<const Topology>(pipeline::CTX_TOPOLOGY_KEY, topo);
      return ComputationResult(true);
    } catch (const std::exception& e) {
      return ComputationResult(ComputationError(std::string("BuildTopology failed: ") + e.what()));
    } catch (...) {
      return ComputationResult(ComputationError("BuildTopology failed")); 
    }
  }
};

} // namespace lahuta::analysis::topology

#endif // LAHUTA_ANALYSIS_TOPOLOGY_KERNEL_HPP
