#ifndef LAHUTA_TASKS_MODEL_WRITE_TASK_HPP
#define LAHUTA_TASKS_MODEL_WRITE_TASK_HPP

#include "gemmi/gz.hpp"
#include "logging.hpp"
#include "mmap/MemoryMapped.h"
#include "models/topology.hpp"
#include <string>
#include <string_view>

namespace lahuta::tasks {

class ModelWriteTask {
public:
  struct result_type {
    bool success;
    std::string file_path;
    ModelParserResult data;
  };

  result_type operator()(std::string_view file_path) const {
    result_type result;
    result.file_path = std::string(file_path);
    result.success = false;

    try {
      if (gemmi::iends_with(result.file_path, ".gz")) {
        gemmi::CharArray buffer = gemmi::MaybeGzipped(result.file_path).uncompress_into_buffer();
        result.data = parse_model(buffer.data(), buffer.size());
        if (!mock_build_model_topology(result.data)) {
          Logger::get_logger()->error("Failed to build topology for file: {}", result.file_path);
          return result;
        }
        result.success = true;
        return result;
      }

      MemoryMapped mm(result.file_path);
      if (!mm.isValid()) {
        Logger::get_logger()->critical("Error opening file: {}", result.file_path);
        return result;
      }
      const char *data = reinterpret_cast<const char *>(mm.getData());
      size_t size = static_cast<size_t>(mm.size());
      result.data = parse_model(data, size);
      if (!mock_build_model_topology(result.data)) {
        Logger::get_logger()->error("Failed to build topology for file: {}", result.file_path);
        return result;
      }
      result.success = true;
    } catch (const std::exception &e) {
      Logger::get_logger()->error("Error processing file {}: {}", result.file_path, e.what());
    }

    return result;
  }
};

} // namespace lahuta::tasks

#endif // LAHUTA_TASKS_MODEL_WRITE_TASK_HPP
