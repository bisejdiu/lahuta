#ifndef LAHUTA_DEFINITIONS_HPP
#define LAHUTA_DEFINITIONS_HPP

#include <array>
#include <string>
#include <unordered_map>
#include <unordered_set>
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

const std::unordered_set<std::string> PositivelyChargedResidues = {"ARG", "HIS", "LYS"};
const std::unordered_set<std::string> NegativelyChargedResidues = {"GLU", "ASP"};

const std::unordered_set<std::string> PolymerNames = {
  "HIS", "ARG", "LYS", "ILE", "PHE", "LEU", "TRP",
  "ALA", "MET", "PRO", "CYS", "ASN", "VAL", "GLY",
  "SER", "GLN", "TYR", "ASP", "GLU", "THR"
};

const std::unordered_set<std::string> BaseNames = {
  "DA", "DC", "DT", "DG", "DI", "DU", "DN"
};

const std::unordered_set<std::string> ProteinBackboneAtoms = {
  "CA",  "C",   "N",   "O",   "O1",  "O2", "OC1",
  "OC2", "OT1", "OT2", "OX1", "OXT", "H",  "H1",
  "H2",  "H3",  "HA",  "HN",  "HXT", "BB"
};

const std::unordered_set<std::string> NucleicBackboneAtoms = {
  "P", "OP1", "OP2", "HOP2", "HOP3", "O2\'", "O3\'",
  "O4\'", "O5\'", "C1\'", "C2\'", "C3\'", "C4\'",
  "C5\'", "H1\'", "H2\'", "H2\'\'", "HO2\'", "H3\'",
  "H4\'", "H5\'", "H5\'\'", "HO3\'", "HO5\'", "O2*",
  "O3*", "O4*", "O5*", "C1*", "C2*", "C3*", "C4*", "C5*"
};

// clang-format on
} // namespace definitions
} // namespace lahuta

#endif // LAHUTA_DEFINITIONS_HPP
