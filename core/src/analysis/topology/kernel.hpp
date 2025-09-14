#pragma once

#include <memory>
#include <string>

#include "analysis/system/records.hpp"
#include "compute/result.hpp"
#include "lahuta.hpp"
#include "pipeline/compute/context.hpp"
#include "pipeline/compute/parameters.hpp"
#include "topology.hpp"
#include <logging.hpp>

// clang-format off
namespace lahuta::analysis::topology {
using namespace pipeline::compute;

struct BuildTopologyKernel {
  static ComputationResult execute(DataContext<PipelineContext, Mut::ReadWrite>& context, const BuildTopologyParams& p) {
    try {
      auto& data = context.data();
      std::shared_ptr<Luni> sys;
      if (data.ctx) sys = std::const_pointer_cast<Luni>(data.ctx->get_object<Luni>("system"));
      if (!sys) return ComputationResult(ComputationError("BuildTopology requires system in context"));

      sys->set_atom_typing_method(p.atom_typing_method);

      // Decide if we're in model fast-path mode
      bool is_model_mode = false;
      if (data.ctx) {
        if (auto mode = data.ctx->get_text("system_mode"); mode && *mode == "model") {
          is_model_mode = true;
        } else {
          // presence of model_data (read from LMDB) implies model mode
          // NOTE: we should not hit this anymore
          auto md = data.ctx->get_object<analysis::system::ModelRecord>("model_data");
          if (md) is_model_mode = true;
        }
      }

      if (is_model_mode) {
        // Ensure topology built via model pathway if not already built
        TopologyBuildingOptions opts{};
        opts.mode = TopologyBuildMode::Model;
        opts.atom_typing_method = p.atom_typing_method;
        (void)sys->build_topology(opts); // idempotent if already built
        sys->enable_only(p.flags);
        auto topo = sys->get_topology_shared();
        if (data.ctx && topo) data.ctx->set_object<Topology>("topology", topo);
        return ComputationResult(true);
      }

      // exec full topology comp
      sys->enable_only(p.flags);

      // edge case: allow disabling all computations while still materializing a Topology object in context.
      if (p.flags == TopologyComputation::None) {
        auto topo0 = sys->get_topology_shared();
        if (!topo0) {
          sys->enable_only(TopologyComputation::None);
          topo0 = sys->get_topology_shared();
        }
        if (data.ctx && topo0) data.ctx->set_object<Topology>("topology", topo0);
        return ComputationResult(true);
      }

      // normal full build
      TopologyBuildingOptions opts{};
      opts.compute_nonstandard_bonds = has_flag(p.flags, TopologyComputation::NonStandardBonds);
      opts.atom_typing_method = p.atom_typing_method;
      if (!sys->build_topology(opts)) return ComputationResult(ComputationError("BuildTopology failed"));
      sys->enable_only(p.flags);
      auto topo = sys->get_topology_shared();
      if (data.ctx && topo) data.ctx->set_object<Topology>("topology", topo);
      return ComputationResult(true);
    } catch (const std::exception& e) {
      return ComputationResult(ComputationError(std::string("BuildTopology failed: ") + e.what()));
    } catch (...) {
      return ComputationResult(ComputationError("BuildTopology failed")); 
    }
  }
};

} // namespace lahuta::analysis::topology
