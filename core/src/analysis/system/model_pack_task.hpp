#ifndef LAHUTA_ANALYSIS_SYSTEM_MODEL_PACK_TASK_HPP
#define LAHUTA_ANALYSIS_SYSTEM_MODEL_PACK_TASK_HPP

#include <memory>
#include <string>
#include <vector>

#include "analysis/system/model_loader.hpp"
#include "analysis/system/records.hpp"
#include "logging/logging.hpp"
#include "pipeline/task/api.hpp"
#include "serialization/formats.hpp"
#include "serialization/serializer.hpp"

namespace lahuta::analysis {
namespace P = lahuta::pipeline;

// ModelPackTask
// Emission on channel `channel_` containing a binary-serialized ModelRecord
class ModelPackTask : public P::ITask {
public:
  explicit ModelPackTask(std::string channel = "db") : channel_(std::move(channel)) {}

  P::TaskResult run(const std::string &item_path, P::TaskContext & /*ctx*/) override {
    try {
      ModelRecord rec = build_model_record(item_path);
      // Thread-local payload buffer to reduce repeated allocations.
      static thread_local std::string tls_payload;
      serialization::Serializer<fmt::binary, ModelRecord>::serialize_into(rec, tls_payload);
      std::string payload;
      payload.swap(tls_payload);
      return P::TaskResult{true, std::vector<P::Emission>{{channel_, std::move(payload)}}};
    } catch (const std::exception &e) {
      Logger::get_logger()->error("[model-pack] Error processing file {}: {}", item_path, e.what());
      return {false, {}};
    }
  }

private:
  std::string channel_;
};

} // namespace lahuta::analysis

#endif // LAHUTA_ANALYSIS_SYSTEM_MODEL_PACK_TASK_HPP
