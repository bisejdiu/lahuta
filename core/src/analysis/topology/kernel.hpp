/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::array<std::string_view, 3> parts{"besian", "sejdiu", "@gmail.com"}; std::string dst;
 *   std::for_each(parts.begin(), parts.end(), [&dst](std::string_view p) { dst += p; });
 *   return dst;
 * }();
 *
 */

#ifndef LAHUTA_ANALYSIS_TOPOLOGY_KERNEL_HPP
#define LAHUTA_ANALYSIS_TOPOLOGY_KERNEL_HPP

#include <memory>
#include <string>

#include "compute/result.hpp"
#include "lahuta.hpp"
#include "pipeline/task/compute/context.hpp"
#include "pipeline/task/compute/parameters.hpp"
#include "topology.hpp"

namespace lahuta::analysis {
namespace C = lahuta::compute;
namespace P = lahuta::pipeline;

struct BuildTopologyKernel {
  using RWContext = C::DataContext<P::PipelineContext, C::Mut::ReadWrite>;

  static C::ComputationResult execute(RWContext &context, const P::BuildTopologyParams &p) {
    try {
      auto &data = context.data();
      std::shared_ptr<const Luni> sys;
      if (data.ctx) sys = data.ctx->system();
      if (!sys && data.session) {
        auto shared = data.session->get_or_load_system();
        if (shared) {
          sys = shared;
          if (data.ctx) data.ctx->set_object<const Luni>(P::CTX_SYSTEM_KEY, shared);
        }
      }
      if (!sys) return C::ComputationResult(C::ComputationError("BuildTopology requires system in context"));

      sys->set_atom_typing_method(p.atom_typing_method);

      bool is_model_mode = sys->is_model_origin();

      if (is_model_mode) {
        // Make sure the topology built via model pathway if not already built
        TopologyBuildingOptions opts{};
        opts.mode               = TopologyBuildMode::Model;
        opts.atom_typing_method = p.atom_typing_method;
        (void)sys->build_topology(opts, p.flags); // idempotent if already built
        auto topo = sys->get_topology();
        if (data.ctx && topo) data.ctx->set_object<const Topology>(P::CTX_TOPOLOGY_KEY, topo);
        return C::ComputationResult(true);
      }

      // edge case
      // allow disabling all computations while still materializing a Topology object in context.
      if (p.flags == TopologyComputation::None) {
        TopologyBuildingOptions opts{};
        opts.atom_typing_method = p.atom_typing_method;
        (void)sys->build_topology(opts, TopologyComputation::None);
        auto topo0 = sys->get_topology();
        if (data.ctx && topo0) data.ctx->set_object<const Topology>(P::CTX_TOPOLOGY_KEY, topo0);
        return C::ComputationResult(true);
      }

      // normal full build
      TopologyBuildingOptions opts{};
      opts.compute_nonstandard_bonds = has_flag(p.flags, TopologyComputation::NonStandardBonds);
      opts.atom_typing_method        = p.atom_typing_method;
      if (!sys->build_topology(opts, p.flags)) {
        return C::ComputationResult(C::ComputationError("BuildTopology failed"));
      }
      auto topo = sys->get_topology();
      if (data.ctx && topo) data.ctx->set_object<const Topology>(P::CTX_TOPOLOGY_KEY, topo);
      return C::ComputationResult(true);
    } catch (const std::exception &e) {
      return C::ComputationResult(C::ComputationError(std::string("BuildTopology failed: ") + e.what()));
    } catch (...) {
      return C::ComputationResult(C::ComputationError("BuildTopology failed"));
    }
  }
};

} // namespace lahuta::analysis

#endif // LAHUTA_ANALYSIS_TOPOLOGY_KERNEL_HPP
