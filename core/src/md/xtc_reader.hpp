#ifndef LAHUTA_MD_XTC_READER_HPP
#define LAHUTA_MD_XTC_READER_HPP

#include <array>
#include <cstdint>
#include <istream>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <rdkit/Geometry/point.h>

#include "md/xdr.hpp"
#include "span.hpp"

// clang-format off
namespace lahuta::md {

class XtcReader {
public:
  explicit XtcReader(std::string_view filename);
  explicit XtcReader(std::istream     &stream);

  struct FramePacket {
    std::shared_ptr<const std::vector<char>> payload;
    std::size_t atom_count = 0;
    std::int64_t step = 0;
    std::optional<double> timestamp_ps;
  };

  void initialize();
  [[nodiscard]] bool read_next_frame();
  // Decode the next frame directly into the provided buffer
  // Returns false on EOF. On success, updates header/time/precision similar to read_next_frame.
  [[nodiscard]] bool read_next_into(RDGeom::POINT3D_VECT &out);
  [[nodiscard]] bool read_next_packet(FramePacket &packet);
  // Random access helpers for seekable streams.
  // - seek_time: position to the first frame with time >= target (or closest when forward_only=false).
  // - seek_step: position to the frame with exact step number if present. Returns false if not found.
  [[nodiscard]] bool seek_time(float target_time, bool forward_only = true);
  [[nodiscard]] bool seek_step(std::int64_t target_step);

  std::size_t natoms() const { return natoms_; }
  std::int64_t current_step() const { return current_header_.step; }
  float current_time() const { return current_header_.time; }
  const std::vector<RDGeom::Point3D> &coords() const { return coords_; }
  // Borrowed zero-copy view, valid until next read
  span<const RDGeom::Point3D> current_coords() const {
    return span<const RDGeom::Point3D>(coords_.data(), coords_.size());
  }
  std::vector<RDGeom::Point3D> move_current_coords();
  const RDGeom::Point3D &periodic_box() const { return box_; }
  // Full triclinic cell 3x3 in Angstroms, row-major
  std::array<float, 9> box_matrix_angstrom() const {
    std::array<float, 9> out = current_header_.box;
    for (float &v : out) {
      v *= 10.0f; // nm -> Angstrom
    }
    return out;
  }
  float current_precision() const { return current_precision_; }

  // Build and query capabilities
  [[nodiscard]] bool is_seekable() const;
  std::streampos file_size_bytes() const;
  [[nodiscard]] bool build_index() { return ensure_index(); }

  // Config
  struct Caps {
    std::uint32_t max_atoms = 20'000'000u;
    std::uint32_t max_bitstream_bytes_per_atom = 256u; // total cap: natoms * factor + 64
    float max_box_length_angstrom = 1.0e6f;
  };
  void set_strict(bool v) { strict_ = v; }
  void set_caps(const Caps &c) { caps_ = c; }
  const Caps &caps() const { return caps_; }

private:
  struct FrameHeader {
    std::int32_t magic   = 0;
    std::uint32_t natoms = 0;
    std::int64_t step    = 0;
    float time;
    std::array<float, 9> box;
  };
  struct FrameIndexEntry {
    std::streampos offset;
    std::int64_t step;
    float time;
  };

  std::string filename_;
  std::shared_ptr<std::istream> stream_ptr_;
  std::istream *stream_;
  std::unique_ptr<internal::XDRReader> xdr_;

  FrameHeader current_header_{};
  RDGeom::Point3D box_{};
  std::vector<RDGeom::Point3D> coords_;
  std::size_t natoms_ = 0;
  bool has_current_ = false;
  float current_precision_ = 0.0f; // -1.0f for uncompressed frames, > 0 for compressed
  bool strict_ = false;
  Caps caps_{};

  // Reusable buffers for compressed decoding
  std::vector<std::byte> bitbuf_;
  std::vector<int> coordbuf_;

  static constexpr unsigned int min_compressed_system_size = 9;

  [[nodiscard]] bool read_frame_header(FrameHeader &header);
  [[nodiscard]] bool read_uncompressed_coords_into(RDGeom::POINT3D_VECT &dest, std::uint32_t atom_count);
  [[nodiscard]] bool read_compressed_coords_into  (RDGeom::POINT3D_VECT &dest, std::uint32_t atom_count);
  [[nodiscard]] bool begin_frame(FrameHeader &header, FramePacket *packet_meta, std::uint32_t &payload_atoms);
  void update_current_index();

  // Seeking utilities
  // Scan forward from current position for the next frame header. Returns header file position or -1.
  [[nodiscard]] std::optional<std::streampos> find_next_frame_start();
  // Check if we are at the start of a header. On success, fills out step/time and returns true.
  [[nodiscard]] std::optional<std::pair<std::int64_t, float>> at_header_start();
  [[nodiscard]] bool ensure_index();
  [[nodiscard]] bool seek_to_index(std::size_t idx);

  std::vector<FrameIndexEntry> frame_index_;
  bool frame_index_ready_ = false;
  std::size_t current_index_ = 0;
};

} // namespace lahuta::md

#endif // LAHUTA_MD_XTC_READER_HPP
