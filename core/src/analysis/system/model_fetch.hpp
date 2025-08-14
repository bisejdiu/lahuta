#pragma once

#include "analysis/system/records.hpp"
#include "compute/dependency.hpp"
#include "db/reader.hpp"
#include "pipeline/compute/context.hpp"
#include "pipeline/compute/dynamic_computation.hpp"
#include "pipeline/compute/parameters.hpp"
#include "serialization/formats.hpp"
#include "serialization/serializer.hpp"

// clang-format off
namespace lahuta::analysis::system {
using namespace lahuta::topology::compute;
using namespace lahuta::pipeline::compute;

// Fetches a model record from LMDB
struct ModelFetchKernel {
  static ComputationResult execute(DataContext<PipelineContext, Mut::ReadWrite>& context, const ModelFetchParams& p) {
    try {
      if (!p.db) return ComputationResult(ComputationError("ModelFetch: database handle is null"));
      auto& data = context.data();
      if (!data.ctx) return ComputationResult(ComputationError("ModelFetch: TaskContext is null"));

      // The input item is the key for LMDB mode
      const std::string& key = data.item_path;
      std::string_view raw;

      // Reuse a thread-local reader per worker
      struct TLReader { std::unique_ptr<LMDBReader> rdr; };
      static thread_local TLReader tls;
      if (!tls.rdr) {
        tls.rdr = std::make_unique<LMDBReader>(p.db->get_env(), p.db->get_dbi());
      }

      if (!tls.rdr->fetch(key, raw)) {
        return ComputationResult(ComputationError("ModelFetch: key not found: " + key));
      }

      ModelRecord obj = serialization::Serializer<fmt::binary, ModelRecord>::deserialize(raw.data(), raw.size());
      data.ctx->set_object<ModelRecord>("model_data", std::make_shared<ModelRecord>(std::move(obj)));
      return ComputationResult(true);

    } catch (const std::exception& e) {
      return ComputationResult(ComputationError(std::string("ModelFetch failed: ") + e.what()));
    } catch (...) {
      return ComputationResult(ComputationError("ModelFetch failed"));
    }
  }
};

// Dynamic label computation with no dependencies
class ModelFetchComputation : public DynamicLabelComputation<ModelFetchParams, ModelFetchKernel, ModelFetchComputation> {
public:
  using Base = DynamicLabelComputation<ModelFetchParams, ModelFetchKernel, ModelFetchComputation>;
  using Base::DynamicLabelComputation;
  using dependencies = UnitComputation;
};

} // namespace lahuta::analysis::system
