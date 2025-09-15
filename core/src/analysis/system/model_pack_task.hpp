#ifndef LAHUTA_ANALYSIS_SYSTEM_MODEL_PACK_TASK_HPP
#define LAHUTA_ANALYSIS_SYSTEM_MODEL_PACK_TASK_HPP

#include <memory>
#include <string>
#include <vector>

#include "analysis/system/records.hpp"
#include "gemmi/gz.hpp"
#include "logging.hpp"
#include "mmap/MemoryMapped.h"
#include "models/parser.hpp"
#include "models/topology.hpp"
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
    WriterRes rec{};
    rec.file_path = item_path;
    rec.success = false;

    try {
      if (gemmi::iends_with(rec.file_path, ".gz")) {
        gemmi::CharArray buffer = gemmi::MaybeGzipped(rec.file_path).uncompress_into_buffer();
        rec.data = parse_model(buffer.data(), buffer.size());
      } else {
        MemoryMapped mm(rec.file_path);
        if (!mm.isValid()) {
          Logger::get_logger()->critical("Error opening file: {}", rec.file_path);
          return {false, {}};
        }
        const char* data = reinterpret_cast<const char*>(mm.getData());
        size_t size = static_cast<size_t>(mm.size());
        rec.data = parse_model(data, size);
      }

      if (!mock_build_model_topology(rec.data)) {
        Logger::get_logger()->error("Failed to build topology for file: {}", rec.file_path);
        // Do not abort, still store record with success=false
      } else {
        rec.success = true;
      }

      std::string payload = serialization::Serializer<fmt::binary, WriterRes>::serialize(rec);
      return TaskResult{true, std::vector<Emission>{{channel_, std::move(payload)}}};

    } catch (const std::exception& e) {
      Logger::get_logger()->error("Error processing file {}: {}", rec.file_path, e.what());
      return {false, {}};
    }
  }

private:
  std::string channel_;
};

} // namespace lahuta::analysis::system

#endif // LAHUTA_ANALYSIS_SYSTEM_MODEL_PACK_TASK_HPP
