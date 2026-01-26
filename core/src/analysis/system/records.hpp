#ifndef LAHUTA_ANALYSIS_SYSTEM_RECORDS_HPP
#define LAHUTA_ANALYSIS_SYSTEM_RECORDS_HPP

#include <string>

#include "models/parser.hpp"

namespace lahuta::analysis {

struct ModelRecord {
  bool success;
  std::string file_path;
  ModelParserResult data;
};

} // namespace lahuta::analysis

#endif // LAHUTA_ANALYSIS_SYSTEM_RECORDS_HPP
