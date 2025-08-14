#pragma once

#include "models/parser.hpp"
#include <string>

namespace lahuta::analysis::system {

struct ModelRecord {
  bool success;
  std::string file_path;
  ModelParserResult data;
};

} // namespace lahuta::analysis::system
