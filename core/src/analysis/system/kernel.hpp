/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::array<std::string_view, 3> parts{"besian", "sejdiu", "@gmail.com"}; std::string s;
 *   for (auto p : parts) s.append(p);
 *   return s;
 * }();
 *
 */

#ifndef LAHUTA_ANALYSIS_SYSTEM_KERNEL_HPP
#define LAHUTA_ANALYSIS_SYSTEM_KERNEL_HPP

#include <memory>
#include <string>

#include "analysis/system/model_loader.hpp"
#include "compute/result.hpp"
#include "lahuta.hpp"
#include "pipeline/task/compute/context.hpp"
#include "pipeline/task/compute/parameters.hpp"
#include "pipeline/task/context.hpp"
#include "pipeline/task/keys.hpp"

namespace lahuta::analysis {
namespace C = lahuta::compute;
namespace P = lahuta::pipeline;

struct SystemReadKernel {
  using RWContext = C::DataContext<P::PipelineContext, C::Mut::ReadWrite>;

  static C::ComputationResult execute(RWContext &context, const P::SystemReadParams &p) {
    try {
      auto &data = context.data();

      std::shared_ptr<const Luni> system;

      if (data.session) {
        system = data.session->get_or_load_system();
      } else if (p.is_model) {
        auto parsed = load_model_parser_result(data.item_path);
        auto s      = Luni::from_model_data(parsed);
        system      = std::make_shared<Luni>(std::move(s));
      } else {
        system = std::make_shared<Luni>(data.item_path);
      }

      if (!system) return C::ComputationResult(C::ComputationError("SystemRead failed: null system"));

      if (data.ctx) {
        data.ctx->set_object<const Luni>(P::CTX_SYSTEM_KEY, system);
      }
      return C::ComputationResult(true);
    } catch (const std::exception &e) {
      return C::ComputationResult(C::ComputationError(std::string("SystemRead failed: ") + e.what()));
    } catch (...) {
      return C::ComputationResult(C::ComputationError("SystemRead failed"));
    }
  }
};

} // namespace lahuta::analysis

#endif // LAHUTA_ANALYSIS_SYSTEM_KERNEL_HPP
