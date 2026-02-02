#ifndef LAHUTA_ANALYSIS_DSSP_KERNEL_HPP
#define LAHUTA_ANALYSIS_DSSP_KERNEL_HPP

// DSSP secondary structure assignment implementation.
//
// This implementation follows the algorithm originally described in:
//   Kabsch W, Sander C. "Dictionary of protein secondary structure: pattern
//   recognition of hydrogen-bonded and geometrical features."
//   Biopolymers. 1983 Dec;22(12):2577-637. doi:10.1002/bip.360221211
//
// Implementation approach informed by PDB-REDO/dssp (https://github.com/PDB-REDO/dssp)
// Copyright (c) 2020 NKI/AVL, Netherlands Cancer Institute
// Licensed under BSD-2-Clause

#include <memory>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "analysis/dssp/hbond.hpp"
#include "analysis/dssp/precompute.hpp"
#include "analysis/dssp/records.hpp"
#include "analysis/dssp/secondary.hpp"
#include "analysis/extract/extract_tasks.hpp"
#include "compute/result.hpp"
#include "pipeline/task/compute/context.hpp"
#include "pipeline/task/compute/parameters.hpp"
#include "pipeline/task/context.hpp"
#include "pipeline/task/emission.hpp"
#include "serialization/serializer.hpp"

namespace lahuta::analysis {
namespace C = lahuta::compute;
namespace P = lahuta::pipeline;

struct DsspKernel {
  static C::ComputationResult execute(C::DataContext<P::PipelineContext, C::Mut::ReadWrite> &context,
                                      const P::DsspParams &p) {
    auto &data = context.data();

    try {
      std::vector<DsspResidue> residues;
      std::string error;

      const P::ModelPayloadSlices *payload = data.ctx ? data.ctx->model_payload().get() : nullptr;
      const Topology *topology             = data.ctx ? data.ctx->topology().get() : nullptr;
      const RDKit::Conformer *conformer    = data.ctx ? data.ctx->conformer().get() : nullptr;

      bool built = false;

      if (conformer && topology) {
        if (build_residues_from_topology(*topology, *conformer, residues, error)) {
          built = true;
        }
      }

      std::string_view sequence;
      if (!built && payload) {
        if (payload->sequence && !payload->sequence->empty()) {
          sequence = *payload->sequence;
        } else if (payload->sequence_view && !payload->sequence_view->data.empty()) {
          sequence = payload->sequence_view->data;
        }

        if (!sequence.empty()) {
          if (payload->positions && !payload->positions->empty()) {
            if (build_residues_from_payload(sequence,
                                            span<const RDGeom::Point3D>(*payload->positions),
                                            residues,
                                            error)) {
              built = true;
            }
          } else if (payload->positions_view && !payload->positions_view->data.empty()) {
            if (build_residues_from_payload(sequence, payload->positions_view->data, residues, error)) {
              built = true;
            }
          }
        }
      }

      if (!built && topology) {
        const auto &conf = topology->conformer();
        if (build_residues_from_topology(*topology, conf, residues, error)) {
          built = true;
        }
      }

      std::shared_ptr<const ModelParserResult> parsed;
      if (!built && data.ctx) {
        parsed = get_cached_model_parser_result(*data.ctx);
      }
      if (!built && parsed && !parsed->sequence.empty() && parsed->coords_size() > 0) {
        sequence = parsed->sequence;
        if (build_residues_from_payload(sequence,
                                        span<const RDGeom::Point3D>(parsed->coords),
                                        residues,
                                        error)) {
          built = true;
        }
      }
      if (!built && error.empty()) {
        if (parsed) {
          error = "cached model parse missing sequence/coords";
        } else {
          error = "no cached model parse; run parse_model or provide topology/payload";
        }
      }

      if (!built) {
        const std::string message = "DSSP residue build failed for '" + data.item_path + "': " + error;
        return C::ComputationResult(C::ComputationError(message));
      }

      compute_backbone_geometry(residues);
      auto ca_pairs = collect_ca_pairs(residues);
      compute_hbond_energies(residues, ca_pairs);
      assign_secondary_structure(residues, ca_pairs, p.pp_stretch_length, p.prefer_pi_helices);

      DsspRecord record;
      record.model_path = data.item_path;
      record.assignments.reserve(residues.size());
      for (const auto &res : residues) {
        record.assignments.push_back(res.assignment);
      }

      auto payload_out = serialization::Serializer<fmt::json, DsspRecord>::serialize(record);

      P::EmissionList out;
      out.push_back(P::Emission{p.channel, std::move(payload_out)});
      return C::ComputationResult(std::move(out));
    } catch (const std::exception &e) {
      return C::ComputationResult(C::ComputationError(std::string("DSSP failed: ") + e.what()));
    } catch (...) {
      return C::ComputationResult(C::ComputationError("DSSP failed"));
    }
  }
};

} // namespace lahuta::analysis

#endif // LAHUTA_ANALYSIS_DSSP_KERNEL_HPP
