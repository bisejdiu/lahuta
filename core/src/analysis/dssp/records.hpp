#ifndef LAHUTA_ANALYSIS_DSSP_RECORDS_HPP
#define LAHUTA_ANALYSIS_DSSP_RECORDS_HPP

#include <string>
#include <string_view>
#include <vector>

#include "models/dssp.hpp"

namespace lahuta::analysis {

inline constexpr std::string_view DsspOutputChannel = "per_residue_dssp";

struct DsspRecord {
  std::string model_path;
  std::vector<DSSPAssignment> assignments;
};

} // namespace lahuta::analysis

#endif // LAHUTA_ANALYSIS_DSSP_RECORDS_HPP
