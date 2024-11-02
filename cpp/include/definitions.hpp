#ifndef LAHUTA_DEFINITIONS_HPP
#define LAHUTA_DEFINITIONS_HPP

#include <array>
#include <string>
#include <unordered_map>
#include <vector>

namespace lahuta {
namespace definitions {
// clang-format off

inline const std::array<std::string, 11>
WaterResidues = {
  "HOH", "W", "SOL", "TIP3",
  "SPC", "H2O", "TIP4", "TIP",
  "DOD", "D3O", "WAT"
};

inline const std::vector<std::string> 
HistidineResidues = {
  "HIS", "HID", "HIE", "HIP"
};

inline const std::unordered_map<std::string, std::vector<int>> 
AromaticResidues = {
  {"PHE", {6}},
  {"TYR", {6}},
  {"HIS", {5}},
  {"TRP", {5, 6}}
};

// clang-format on
} // namespace definitions
} // namespace lahuta

#endif // LAHUTA_DEFINITIONS_HPP
