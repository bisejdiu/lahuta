#ifndef LAHUTA_PIPELINE_DATA_FRAME_HPP
#define LAHUTA_PIPELINE_DATA_FRAME_HPP

#include <cstddef>
#include <memory>
#include <optional>
#include <string>

#include <rdkit/Geometry/point.h>

namespace lahuta::pipeline {

struct CoordinatesView {
  RDGeom::POINT3D_VECT positions;
  std::shared_ptr<const RDGeom::POINT3D_VECT> shared_positions;

  const RDGeom::POINT3D_VECT &data() const { return shared_positions ? *shared_positions : positions; }
};

struct FrameMetadata {
  std::string session_id;
  std::uint64_t conformer_id = 0;
  std::optional<double> timestamp_ps;
  std::optional<std::string> source_file; // trajectory file for MD data
};

struct FrameHandle {
  virtual ~FrameHandle() = default;

  virtual std::size_t index() const noexcept                  = 0;
  virtual std::optional<double> timestamp_ps() const noexcept = 0;
  virtual CoordinatesView load_coordinates()                  = 0;
};

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_DATA_FRAME_HPP
