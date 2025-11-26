#ifndef LAHUTA_DB_MODEL_PAYLOAD_HPP
#define LAHUTA_DB_MODEL_PAYLOAD_HPP

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string_view>
#include <type_traits>

// clang-format off
namespace lahuta {

enum class ModelPayloadField : std::uint8_t {
  Sequence = 0,
  TaxonomyId,
  OrganismScientific,
  Coordinates,
  PlddtCategories,
  DsspAssignments,
  Count
};

struct ModelFieldSlice {
  std::uint32_t offset = 0;
  std::uint32_t length = 0;
};

struct ModelPayloadHeader {
  std::uint32_t magic    = 0;
  std::uint16_t version  = 0;
  std::uint8_t  flags    = 0;
  std::uint8_t  reserved = 0;
  std::uint32_t sequence_length = 0;
  std::uint32_t coord_count = 0;
  ModelFieldSlice slices[static_cast<std::size_t>(ModelPayloadField::Count)]{};
};

inline constexpr std::uint32_t ModelPayloadMagic = 0x4C4D444D; // LMDM
inline constexpr std::uint16_t ModelPayloadVersion = 1;
inline constexpr std::uint8_t  ModelPayloadFlagCompressed = 0x01;
inline constexpr std::size_t   CoordinateStride = 3 * sizeof(float);

class ModelPayloadView {
public:
  explicit ModelPayloadView(std::string_view blob) : blob_(blob) {
    if (blob_.size() < sizeof(ModelPayloadHeader)) { throw std::runtime_error("Model payload too small"); }
    std::memcpy(&header_, blob_.data(), sizeof(header_));

    if (header_.magic   != ModelPayloadMagic)       { throw std::runtime_error("Model payload missing magic header"); }
    if (header_.version != ModelPayloadVersion)     { throw std::runtime_error("Unsupported model payload version"); }
    if (header_.flags & ModelPayloadFlagCompressed) { throw std::runtime_error("Compressed model payloads are not supported yet"); }

    for (const auto& slice : header_.slices) {
      validate_slice(slice);
    }
  }

  std::uint32_t coord_count() const noexcept { return header_.coord_count; }
  std::uint32_t sequence_length() const noexcept { return header_.sequence_length; }
  ModelFieldSlice coord_slice() const noexcept {
    return header_.slices[static_cast<std::size_t>(ModelPayloadField::Coordinates)];
  }
  std::string_view blob() const noexcept { return blob_; }

  std::string_view sequence()            const { return field(ModelPayloadField::Sequence); }
  std::string_view taxonomy_id()         const { return field(ModelPayloadField::TaxonomyId); }
  std::string_view organism_scientific() const { return field(ModelPayloadField::OrganismScientific); }

  std::string_view coords_bytes() const {
    const auto slice = coord_slice();
    if (slice.length == 0) return {};
    validate_slice(slice);
    const std::size_t payload_bytes = static_cast<std::size_t>(header_.coord_count) * CoordinateStride;
    if (payload_bytes == 0) {
      throw std::runtime_error("Model payload coordinate slice present but coord_count is zero");
    }
    if (slice.length < payload_bytes) {
      throw std::runtime_error("Model payload coordinate slice too small");
    }
    return std::string_view(blob_.data() + slice.offset, payload_bytes);
  }
  std::string_view plddt_bytes() const { return field(ModelPayloadField::PlddtCategories); }
  std::string_view dssp_bytes()  const { return field(ModelPayloadField::DsspAssignments); }

private:
  std::string_view field(ModelPayloadField f) const {
    const auto& slice = header_.slices[static_cast<std::size_t>(f)];
    if (slice.length == 0) return {};
    validate_slice(slice);
    return std::string_view(blob_.data() + slice.offset, slice.length);
  }

  void validate_slice(const ModelFieldSlice& slice) const {
    if (slice.length == 0) return;
    const std::size_t begin = static_cast<std::size_t>(slice.offset);
    const std::size_t end = begin + static_cast<std::size_t>(slice.length);
    if (end > blob_.size()) {
      throw std::runtime_error("Model payload slice exceeds buffer");
    }
    if (begin < sizeof(ModelPayloadHeader)) {
      throw std::runtime_error("Model payload slice overlaps header");
    }
  }

  ModelPayloadHeader header_{};
  std::string_view blob_;
};

static_assert(std::is_trivially_copyable_v<ModelPayloadHeader>, "Header must be trivially copyable");

} // namespace lahuta

#endif // LAHUTA_DB_MODEL_PAYLOAD_HPP
