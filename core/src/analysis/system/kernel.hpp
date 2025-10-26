#ifndef LAHUTA_ANALYSIS_SYSTEM_KERNEL_HPP
#define LAHUTA_ANALYSIS_SYSTEM_KERNEL_HPP

#include <memory>
#include <string>

#include "analysis/system/model_loader.hpp"
#include "compute/result.hpp"
#include "lahuta.hpp"
#include "models/metadata.hpp"
#include "pipeline/compute/context.hpp"
#include "pipeline/compute/parameters.hpp"
#include "pipeline/dynamic/keys.hpp"

// clang-format off
namespace lahuta::analysis::system {
using namespace lahuta::pipeline::compute;

inline void publish_model_metadata(pipeline::dynamic::TaskContext* ctx, const ModelMetadata& meta) {
  if (!ctx || meta.empty()) return;
  ctx->set_object<ModelMetadata>(pipeline::CTX_MODEL_METADATA_KEY, std::make_shared<ModelMetadata>(meta));
}

struct SystemReadKernel {
  static ComputationResult execute(DataContext<PipelineContext, Mut::ReadWrite>& context, const SystemReadParams& p) {
    try {
      auto& data = context.data();

      std::shared_ptr<const Luni> sys;

      if (data.session) {
        sys = data.session->get_or_load_system();
        if (auto meta = data.session->model_metadata()) {
          publish_model_metadata(data.ctx, *meta);
        }
      } else if (p.is_model) {
        auto parsed = load_model_parser_result(data.item_path);
        publish_model_metadata(data.ctx, parsed.metadata);
        auto s = Luni::from_model_data(parsed);
        sys = std::make_shared<Luni>(std::move(s));
      } else {
        sys = std::make_shared<Luni>(data.item_path);
      }

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
