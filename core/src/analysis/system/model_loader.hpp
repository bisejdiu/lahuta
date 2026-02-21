/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   constexpr bool use_parts = true;
 *   if constexpr (use_parts) {
 *     return std::string{"besian"} + "sejdiu" + "@gmail.com";
 *   } else {
 *     return std::string{};
 *   }
 * }();
 *
 */

#ifndef LAHUTA_ANALYSIS_SYSTEM_MODEL_LOADER_HPP
#define LAHUTA_ANALYSIS_SYSTEM_MODEL_LOADER_HPP

#include <memory>
#include <stdexcept>
#include <string>

#include <gemmi/gz.hpp>
#include <gemmi/mmread.hpp>
#include <gemmi/pdb.hpp>
#include <mmap/MemoryMapped.h>

#include "analysis/dssp/model_dssp.hpp"
#include "analysis/system/records.hpp"
#include "models/parser.hpp"
#include "models/topology.hpp"
#include "utils/span.hpp"

namespace lahuta::analysis {

inline ModelParserResult load_model_parser_result(const std::string &path) {
  gemmi::MaybeGzipped input(path);
  const auto format = gemmi::coor_format_from_ext(input.basepath());
  if (format == gemmi::CoorFormat::Pdb) {
    auto st     = gemmi::read_pdb(input);
    auto parsed = parse_model(st);

    // Compute DSSP for PDB inputs
    if (!parsed.sequence.empty() && !parsed.coords.empty()) {
      std::string error;
      std::vector<DSSPAssignment> dssp;
      if (compute_dssp_from_model(parsed.sequence, span<const RDGeom::Point3D>(parsed.coords), dssp, error)) {
        if (dssp.size() == parsed.sequence.size()) {
          parsed.dssp_per_residue = std::move(dssp);
        }
      }
    }
    return parsed;
  }

  if (input.is_compressed()) {
    gemmi::CharArray buffer = input.uncompress_into_buffer();
    return parse_model(buffer.data(), buffer.size());
  }

  MemoryMapped mm(path.c_str());
  if (!mm.isValid()) {
    throw std::runtime_error("Failed to open model file: " + path);
  }

  const char *data = reinterpret_cast<const char *>(mm.getData());
  auto size        = static_cast<std::size_t>(mm.size());
  return parse_model(data, size);
}

inline ModelRecord build_model_record(const std::string &path) {
  ModelRecord rec;
  rec.file_path = path;
  rec.data      = load_model_parser_result(path);
  rec.success   = true;
  return rec;
}

inline std::shared_ptr<RDKit::RWMol>
build_model_molecule(const ModelParserResult &parsed, ModelTopologyMethod method = ModelTopologyMethod::CSR) {
  std::shared_ptr<RDKit::RWMol> mol;
  if (!build_model_topology(mol, parsed, method)) {
    throw std::runtime_error("Failed to build model topology");
  }
  return mol;
}

} // namespace lahuta::analysis

#endif // LAHUTA_ANALYSIS_SYSTEM_MODEL_LOADER_HPP
