#ifndef LAHUTA_PIPELINE_FRAME_HPP
#define LAHUTA_PIPELINE_FRAME_HPP

#include <cstddef>
#include <memory>
#include <optional>
#include <string>

#include <rdkit/Geometry/point.h>

namespace lahuta {

struct CoordinatesView {
  RDGeom::POINT3D_VECT positions;
  std::shared_ptr<const RDGeom::POINT3D_VECT> shared_positions;

  const RDGeom::POINT3D_VECT &data() const { return shared_positions ? *shared_positions : positions; }
};

struct FrameMetadata {
  std::string session_id;
  std::uint64_t conformer_id = 0;
  std::optional<double> timestamp_ps;
};

struct FrameHandle {
  virtual ~FrameHandle() = default;
  virtual std::size_t index() const noexcept = 0;
  virtual std::optional<double> timestamp_ps() const noexcept = 0;
  virtual CoordinatesView load_coordinates() = 0;
};

} // namespace lahuta

#endif // LAHUTA_PIPELINE_FRAME_HPP
