#pragma once

#include <istream>
#include <memory>
#include <string_view>

#include <GraphMol/RWMol.h>

namespace lahuta::md {

// Read a GRO structure file and return an RDKit::RWMol with a single conformer.
// - Coordinates are converted from nm to A.
// - Atom elements are inferred from atom/residue names using a heuristic mapping.
// - Box line is consumed but ignored (no periodicity information preserved).
// - Velocities are not parsed (coordinates only).
// - Single-frame only (no multi-frame trajectory support).
// - Title line is read but time/step parsing is omitted.
// - Integer field parsing failures default to 0 (resid, atom_id).
std::shared_ptr<RDKit::RWMol> read_gro_to_rwmol(std::string_view filename);

// Internal helper primarily for testing.
std::shared_ptr<RDKit::RWMol> read_gro_stream_to_rwmol(std::istream &stream);

} // namespace lahuta::md
