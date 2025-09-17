#pragma once

#include <string>

#include "models/parser.hpp"

namespace lahuta::analysis::system {

struct ModelRecord {
  bool success;
  std::string file_path;
  ModelParserResult data;
};

} // namespace lahuta::analysis::system
