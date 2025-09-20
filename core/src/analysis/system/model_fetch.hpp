#pragma once

#include <memory>

#include "analysis/system/records.hpp"
#include "db/reader.hpp"
#include "pipeline/compute/context.hpp"
#include "pipeline/compute/parameters.hpp"
#include "serialization/formats.hpp"
#include "serialization/serializer.hpp"

// clang-format off
namespace lahuta::analysis::system {
using namespace lahuta::pipeline::compute;

// Fetches a model record from LMDB
struct ModelFetchKernel {
  static ComputationResult execute(DataContext<PipelineContext, Mut::ReadWrite>& context, const ModelFetchParams& p) {
    try {
      if (!p.db) return ComputationResult(ComputationError("ModelFetch: database handle is null"));
      auto& data = context.data();
      if (!data.ctx) return ComputationResult(ComputationError("ModelFetch: TaskContext is null"));

      const std::string& key = data.item_path;
      std::string_view raw;

      // Reuse a thread-local reader per worker, rebinding when the database changes.
      // `bound_db` keeps a weak reference to the LMDBDatabase currently backing the reader so
      //   we can detect when the shared handle expires or a different database is provided.
      // `rdr` owns the actual LMDBReader instance that is reused until bound_db no longer
      //   matches the requested database, at which point it is recreated against the new env/DBI.
      struct TLReader {
        std::weak_ptr<LMDBDatabase> bound_db;
        std::unique_ptr<LMDBReader> rdr;
      };
      static thread_local TLReader tls;

      auto db = p.db;
      auto locked = tls.bound_db.lock();
      if (!tls.rdr || locked.get() != db.get()) {
        tls.rdr = std::make_unique<LMDBReader>(db->get_env(), db->get_dbi());
        tls.bound_db = db;
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

} // namespace lahuta::analysis::system
