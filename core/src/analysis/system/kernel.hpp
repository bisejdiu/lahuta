#ifndef LAHUTA_ANALYSIS_SYSTEM_KERNEL_HPP
#define LAHUTA_ANALYSIS_SYSTEM_KERNEL_HPP

#include <memory>
#include <string>

#include "compute/result.hpp"
#include "lahuta.hpp"
#include "pipeline/compute/context.hpp"
#include "pipeline/compute/parameters.hpp"

// clang-format off
namespace lahuta::analysis::system {
using namespace lahuta::pipeline::compute;

struct SystemReadKernel {
  static ComputationResult execute(DataContext<PipelineContext, Mut::ReadWrite>& context, const SystemReadParams& p) {
    try {
      auto& data = context.data();

      auto sys = [&data, &p]() -> std::shared_ptr<const Luni> {

        if (data.session) return data.session->get_or_load_system();
        if (!p.is_model)  return std::make_shared<Luni>(data.item_path);

        auto s = Luni::from_model_file(data.item_path);
        return std::make_shared<Luni>(std::move(s));
      }();

      if (!sys) return ComputationResult(ComputationError("SystemRead failed: null system"));

      if (data.ctx) {
        data.ctx->set_object<const Luni>(pipeline::CTX_SYSTEM_KEY, sys);
      }
      return ComputationResult(true);
    } catch (const std::exception& e) {
      return ComputationResult(ComputationError(std::string("SystemRead failed: ") + e.what()));
    } catch (...) {
      return ComputationResult(ComputationError("SystemRead failed"));
    }
  }
};

} // namespace lahuta::analysis::system

#endif // LAHUTA_ANALYSIS_SYSTEM_KERNEL_HPP
