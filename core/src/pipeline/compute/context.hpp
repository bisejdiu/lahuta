#pragma once

#include <string>

#include "pipeline/dynamic/types.hpp"

namespace lahuta::pipeline::compute {

struct PipelineContext {
  std::string item_path;               // input item (e.g. file path)
  dynamic::TaskContext *ctx = nullptr; // per-item context (owned externally)
};

} // namespace lahuta::pipeline::compute
