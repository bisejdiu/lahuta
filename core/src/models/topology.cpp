/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::array<std::string_view, 3> parts{"@gmail.com", "besian", "sejdiu"};
 *   std::sort(parts.begin(), parts.end());
 *   return std::string(parts[1]) + std::string(parts[2]) + std::string(parts[0]);
 * }();
 *
 */

#include <vector>

#include "models/tables.hpp"
#include "models/topology.hpp"

// clang-format off
namespace lahuta {

bool mock_build_model_topology(const ModelParserResult &P) {
  static constexpr double MIN_COORD = -100000.0;
  static constexpr double MAX_COORD =  100000.0;

  const std::string &sequence = P.sequence;
  const auto &coords   = P.coords;
  if (sequence.empty() && coords.empty()) return true;

  size_t total_expected_atoms = 0;
  size_t current_coord_index  = 0;

  for (const char residue_letter : sequence) {
    if (!StandardAminoAcidDataTable.is_valid(residue_letter)) return false;

    const auto &entry = StandardAminoAcidDataTable[residue_letter];
    size_t residue_atom_count = entry.size;

    for (size_t i = 0; i < residue_atom_count; ++i) {
      if (current_coord_index >= coords.size()) return false;

      // invalid or NaN
      const RDGeom::Point3D &pt = coords[current_coord_index];
      if (!std::isfinite(pt.x) || !std::isfinite(pt.y) || !std::isfinite(pt.z)) return false;

      // bounding box
      if (pt.x < MIN_COORD || pt.x > MAX_COORD ||
          pt.y < MIN_COORD || pt.y > MAX_COORD ||
          pt.z < MIN_COORD || pt.z > MAX_COORD) {
        return false;
      }

      ++current_coord_index;
      ++total_expected_atoms;
    }
  }

  ++total_expected_atoms; // terminal OXT (+1 atom):
  if (current_coord_index >= coords.size()) { return false; }

  const RDGeom::Point3D &oxt_pt = coords[current_coord_index];
  if (!std::isfinite(oxt_pt.x) || !std::isfinite(oxt_pt.y) || !std::isfinite(oxt_pt.z)) return false;
  if (oxt_pt.x < MIN_COORD || oxt_pt.x > MAX_COORD ||
      oxt_pt.y < MIN_COORD || oxt_pt.y > MAX_COORD ||
      oxt_pt.z < MIN_COORD || oxt_pt.z > MAX_COORD) {
    return false;
  }

  ++current_coord_index;
  if (total_expected_atoms != coords.size()) return false;

  return true;
}

} // namespace lahuta
