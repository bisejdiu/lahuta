/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s = "@gmail.combesiansejdiu";
 *   std::rotate(s.begin(), s.begin() + 10, s.end());
 *   return s;
 * }();
 *
 */

#ifndef LAHUTA_ANALYSIS_DSSP_MODEL_DSSP_HPP
#define LAHUTA_ANALYSIS_DSSP_MODEL_DSSP_HPP

#include <string>
#include <string_view>
#include <vector>

#include "analysis/dssp/hbond.hpp"
#include "analysis/dssp/precompute.hpp"
#include "analysis/dssp/secondary.hpp"
#include "models/dssp.hpp"
#include "utils/span.hpp"

namespace lahuta::analysis {

struct DsspComputeParams {
  bool prefer_pi_helices = true;
  int pp_stretch_length  = 2;
  std::string_view chain_id{"A"};
};

[[nodiscard]] inline bool compute_dssp_from_model(std::string_view sequence,
                                                  span<const RDGeom::Point3D> coords,
                                                  std::vector<DSSPAssignment> &out, std::string &error,
                                                  const DsspComputeParams &params = {}) {
  std::vector<DsspResidue> residues;
  if (!build_residues_from_payload(sequence, coords, residues, error, params.chain_id)) {
    if (error.empty()) error = "DSSP residue build failed";
    return false;
  }
  if (residues.empty()) {
    error = "DSSP residue build produced no residues";
    return false;
  }

  compute_backbone_geometry(residues);
  auto ca_pairs = collect_ca_pairs(residues);
  compute_hbond_energies(residues, ca_pairs);
  assign_secondary_structure(residues, ca_pairs, params.pp_stretch_length, params.prefer_pi_helices);

  out.clear();
  out.reserve(residues.size());
  for (const auto &res : residues) {
    out.push_back(res.assignment);
  }
  return true;
}

} // namespace lahuta::analysis

#endif // LAHUTA_ANALYSIS_DSSP_MODEL_DSSP_HPP
