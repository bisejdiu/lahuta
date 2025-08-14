#pragma once

#include <memory>
#include <string>

#include "compute/dependency.hpp"
#include "compute/result.hpp"
#include "pipeline/compute/context.hpp"
#include "pipeline/compute/dynamic_computation.hpp"
#include "pipeline/compute/parameters.hpp"
#include "topology.hpp"
#include "topology/compute.hpp"
#include <analysis/topology/computation.hpp>
#include <logging.hpp>

//
// EnsureTypingKernel: Makes sure that the current Topology's atom typing mode matches a desired
// provider-specific mode before downstream analyses run.
//
// - desired:  the requested atom typing mode (ContactComputerType::{Molstar, Arpeggio}).
//             When desired==None, the kernel is a no-op.
// - sentinel: a TaskContext text key ("atom_typing_mode") used as a per-item
//             first-touch marker. The first EnsureTyping that runs sets it to the
//             chosen mode. Subsequent EnsureTyping kernels compare against it to
//             avoid flip-flopping when multiple providers are present.
//
// Practical behavior with contact providers:
// - MolStar contacts: desired=Molstar. If current typing is already MolStar, no-op.
//                     Otherwise, retype to MolStar, set sentinel="molstar".
// - Arpeggio contacts: desired=Arpeggio. If current typing is MolStar, retype to
//                      Arpeggio, set sentinel="arpeggio". If already Arpeggio, no-op.
//
// A correctness guard in ContactsKernel also ensures the right typing at point-of-use.
//
namespace lahuta::analysis::topology {
using namespace lahuta::topology::compute;

// clang-format off
struct EnsureTypingKernel {
  static ComputationResult execute(DataContext<PipelineContext, Mut::ReadWrite>& context, const EnsureTypingParams& p) {
    auto& data = context.data();

    try {
      if (p.desired == ContactComputerType::None) {
        Logger::get_logger()->info("EnsureTyping: desired=None (no-op)");
        return ComputationResult(true);
      }

      if (!data.ctx) return ComputationResult(ComputationError("EnsureTyping requires TaskContext"));

      // Retrieve current topology and compute current typing mode
      auto top_c = data.ctx->get_object<Topology>("topology");
      if (!top_c) return ComputationResult(ComputationError("EnsureTyping requires topology in context"));

      bool use_molstar = true; // default to MolStar when unknown
      {
        using namespace lahuta::topology;
        const auto& label  = AtomTypingComputation<>::label;
        auto& eng = const_cast<Topology&>(*top_c).get_engine();
        auto* params = eng.get_parameters<AtomTypingParams>(label);
        if (params) use_molstar = params->use_molstar;
      }

      const std::string* sentinel = data.ctx->get_text("atom_typing_mode");
      auto desired_label = (p.desired == ContactComputerType::Molstar) ? std::string("molstar") : std::string("arpeggio");

      auto ensure_now = [&](std::shared_ptr<Topology> top_mut) {
        if (p.desired == ContactComputerType::Molstar) {
          Logger::get_logger()->info("EnsureTyping: retyping to molstar");
          top_mut->assign_molstar_typing();
        } else {
          Logger::get_logger()->info("EnsureTyping: retyping to arpeggio");
          top_mut->assign_arpeggio_atom_types();
        }
      };

      // First-touch semantics via sentinel
      if (!sentinel) {
        data.ctx->set_text("atom_typing_mode", desired_label);
        bool desired_is_mol = (p.desired == ContactComputerType::Molstar);
        if (use_molstar != desired_is_mol) {
          auto top_mut = std::const_pointer_cast<Topology>(top_c);
          ensure_now(std::move(top_mut));
        } else {
          Logger::get_logger()->info("EnsureTyping: already in desired mode ({}), no-op", desired_label);
        }
        return ComputationResult(true);
      }

      // If sentinel matches our desired mode, enforce consistency, else no-op
      if (*sentinel == desired_label) {
        bool desired_is_mol = (p.desired == ContactComputerType::Molstar);
        if (use_molstar != desired_is_mol) {
          auto top_mut = std::const_pointer_cast<Topology>(top_c);
          ensure_now(std::move(top_mut));
        } else {
          Logger::get_logger()->info("EnsureTyping: mode '{}' already set and aligned", desired_label);
        }
      } else {
        Logger::get_logger()->info("EnsureTyping: sentinel='{}' differs from desired='{}' - no-op", *sentinel, desired_label);
      }

      return ComputationResult(true);
    } catch (const std::exception& e) {
      return ComputationResult(ComputationError(std::string("EnsureTyping failed: ") + e.what()));
    } catch (...) {
      return ComputationResult(ComputationError("EnsureTyping failed"));
    }
  }
};

class EnsureTypingComputation : public DynamicLabelComputation<EnsureTypingParams, EnsureTypingKernel, EnsureTypingComputation> {
public:
  using Base = DynamicLabelComputation<EnsureTypingParams, EnsureTypingKernel, EnsureTypingComputation>;
  using Base::DynamicLabelComputation;
  using dependencies = Dependencies<Dependency<BuildTopologyComputation, void>>;
};

} // namespace lahuta::analysis::topology
