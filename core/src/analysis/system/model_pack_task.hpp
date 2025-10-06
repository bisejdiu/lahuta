#ifndef LAHUTA_ANALYSIS_SYSTEM_MODEL_PACK_TASK_HPP
#define LAHUTA_ANALYSIS_SYSTEM_MODEL_PACK_TASK_HPP

#include <memory>
#include <string>
#include <vector>

#include "analysis/system/model_loader.hpp"
#include "analysis/system/records.hpp"
#include "logging.hpp"
#include "pipeline/dynamic/types.hpp"
#include "serialization/formats.hpp"
#include "serialization/serializer.hpp"

// clang-format off
namespace lahuta::analysis::system {
using namespace lahuta::pipeline::dynamic;

// ModelPackTask
// Emission on channel `channel_` containing a binary-serialized ModelRecord
class ModelPackTask : public ITask {
public:
  explicit ModelPackTask(std::string channel = "db") : channel_(std::move(channel)) {}

  TaskResult run(const std::string& item_path, TaskContext& /*ctx*/) override {
    using WriterRes = analysis::system::ModelRecord;
    try {
      WriterRes rec = build_model_record(item_path);
      // Thread-local payload buffer to reduce repeated allocations.
      static thread_local std::string tls_payload;
      serialization::Serializer<fmt::binary, WriterRes>::serialize_into(rec, tls_payload);
      std::string payload;
      payload.swap(tls_payload);
      return TaskResult{true, std::vector<Emission>{{channel_, std::move(payload)}}};
    } catch (const std::exception& e) {
      Logger::get_logger()->error("Error processing file {}: {}", item_path, e.what());
      return {false, {}};
    }
  }

private:
  std::string channel_;
};

} // namespace lahuta::analysis::system

#endif // LAHUTA_ANALYSIS_SYSTEM_MODEL_PACK_TASK_HPP
