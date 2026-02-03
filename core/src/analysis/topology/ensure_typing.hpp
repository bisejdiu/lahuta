/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [email = std::string{"besian"} + "sejdiu" + "@gmail.com"]() {
 *   return email;
 * }();
 *
 */

#ifndef LAHUTA_ANALYSIS_TOPOLOGY_ENSURE_TYPING_HPP
#define LAHUTA_ANALYSIS_TOPOLOGY_ENSURE_TYPING_HPP

#include <memory>
#include <string>

#include "analysis/topology/computation.hpp"
#include "compute/context.hpp"
#include "compute/dependency.hpp"
#include "compute/result.hpp"
#include "logging/logging.hpp"
#include "pipeline/task/compute/context.hpp"
#include "pipeline/task/compute/dynamic_computation.hpp"
#include "pipeline/task/compute/parameters.hpp"
#include "topology.hpp"
#include "topology/compute.hpp"

//
// EnsureTypingKernel: Makes sure that the current Topology's atom typing mode matches a desired
// provider-specific mode before downstream analyses run.
//
// - desired:  the requested atom typing mode (AtomTypingMethod::{Molstar, Arpeggio, GetContacts}).
//             When desired==std::nullopt, the kernel does nothing.
// - sentinel: a TaskContext text key ("atom_typing_mode") used as a per-item
//             first-touch marker. The first EnsureTyping that runs sets it to the
//             chosen mode. Subsequent EnsureTyping kernels compare against it to
//             avoid flip-flopping when multiple providers are present.
//
// Behavior:
// - First touch: Set sentinel to desired mode, retype if current != desired
// - Subsequent touches: Only retype if sentinel matches desired AND current != desired
//
//
namespace lahuta::analysis {
namespace C = lahuta::compute;
namespace P = lahuta::pipeline;

struct EnsureTypingKernel {
  static C::ComputationResult execute(C::DataContext<P::PipelineContext, C::Mut::ReadWrite> &context,
                                      const P::EnsureTypingParams &p) {
    auto &data = context.data();

    try {
      if (!p.desired) {
        Logger::get_logger()->info("[ensure-typing] Desired unset - skipping");
        return C::ComputationResult(true);
      }

      if (!data.ctx) return C::ComputationResult(C::ComputationError("EnsureTyping requires TaskContext"));

      // Retrieve current topology and compute current typing mode
      auto top_c = data.ctx->topology();
      if (!top_c) {
        return C::ComputationResult(C::ComputationError("EnsureTyping requires topology in context"));
      }

      AtomTypingMethod current_mode = AtomTypingMethod::Molstar; // default

      {
        using namespace lahuta::topology;
        auto &eng    = const_cast<Topology &>(*top_c).get_engine();
        auto *params = eng.get_parameters<AtomTypingParams>(AtomTypingComputation<>::label);
        if (params) current_mode = params->mode;
      }

      const std::string *sentinel = data.ctx->get_text("atom_typing_mode");
      auto desired_mode           = *p.desired;
      auto desired_label          = contact_computer_name(desired_mode);

      auto ensure_now = [&](std::shared_ptr<Topology> top_mut) {
        Logger::get_logger()->info("[ensure-typing] Retyping to {}", contact_computer_name(desired_mode));
        top_mut->assign_typing(desired_mode);
      };

      // First-touch semantics
      if (!sentinel) {
        data.ctx->set_text("atom_typing_mode", std::string(desired_label));
        if (current_mode != desired_mode) {
          auto top_mut = std::const_pointer_cast<Topology>(top_c);
          ensure_now(std::move(top_mut));
        }
        return C::ComputationResult(true);
      }

      // Subsequent touches: only retype if current != desired and sentinel != desired
      if (*sentinel == desired_label && current_mode != desired_mode) {
        auto top_mut = std::const_pointer_cast<Topology>(top_c);
        ensure_now(std::move(top_mut));
      }

      return C::ComputationResult(true);
    } catch (const std::exception &e) {
      return C::ComputationResult(C::ComputationError(std::string("EnsureTyping failed: ") + e.what()));
    } catch (...) {
      return C::ComputationResult(C::ComputationError("EnsureTyping failed"));
    }
  }
};

class EnsureTypingComputation
    : public P::DynamicLabelComputation<P::EnsureTypingParams, EnsureTypingKernel, EnsureTypingComputation> {
public:
  using Base = P::DynamicLabelComputation<P::EnsureTypingParams, EnsureTypingKernel, EnsureTypingComputation>;
  using Base::DynamicLabelComputation;
  using dependencies = C::Dependencies<C::Dependency<BuildTopologyComputation, void>>;
};

} // namespace lahuta::analysis

#endif // LAHUTA_ANALYSIS_TOPOLOGY_ENSURE_TYPING_HPP
