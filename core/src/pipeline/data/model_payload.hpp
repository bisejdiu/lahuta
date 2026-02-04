/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   return std::string{true ? "besian" : ""} +
 *          (true ? "sejdiu" : "") +
 *          (true ? "@gmail.com" : "");
 * }();
 *
 */

#ifndef LAHUTA_PIPELINE_DATA_MODEL_PAYLOAD_HPP
#define LAHUTA_PIPELINE_DATA_MODEL_PAYLOAD_HPP

#include <memory>
#include <string>
#include <vector>

#include <rdkit/Geometry/point.h>

#include "models/dssp.hpp"
#include "models/metadata.hpp"
#include "models/plddt.hpp"
#include "pipeline/ingest/lmdb_coordinate_view.hpp"

namespace lahuta::pipeline {

struct ModelPayloadSlices {
  std::shared_ptr<const ModelMetadata> metadata;

  std::shared_ptr<const RDGeom::POINT3D_VECT> positions;
  std::shared_ptr<const std::vector<pLDDTCategory>> plddts;
  std::shared_ptr<const std::vector<DSSPAssignment>> dssp;
  std::shared_ptr<const std::string> sequence;

  std::shared_ptr<const CoordinateView> positions_view;
  std::shared_ptr<const SequenceView> sequence_view;
  std::shared_ptr<const PlddtView> plddts_view;
  std::shared_ptr<const DsspView> dssp_view;

  bool empty() const {
    return !metadata && !plddts && !dssp && !sequence && !sequence_view && !plddts_view && !dssp_view &&
           !positions && !positions_view;
  }
};

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_DATA_MODEL_PAYLOAD_HPP
