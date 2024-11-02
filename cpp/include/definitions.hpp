#ifndef LAHUTA_DEFINITIONS_HPP
#define LAHUTA_DEFINITIONS_HPP

#include <array>
#include <string>
#include <vector>

namespace lahuta {

namespace definitions {

// Residue names that represent water molecules
inline const std::array<std::string, 11> WaterResidueNames = {
    "HOH", "W", "SOL", "TIP3", "SPC", "H2O", "TIP4", "TIP", "DOD", "D3O", "WAT"};

inline const std::vector<std::string> HistidineResNames = {"HIS", "HID", "HIE", "HIP"};

} // namespace definitions
} // namespace lahuta

#endif // LAHUTA_DEFINITIONS_HPP
