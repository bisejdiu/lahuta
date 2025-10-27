#ifndef LAHUTA_SERIALIZATION_MODEL_SERIALIZER_HPP
#define LAHUTA_SERIALIZATION_MODEL_SERIALIZER_HPP

#include <cstring>
#include <stdexcept>
#include <string_view>

#include "analysis/system/records.hpp"
#include "db/model_payload.hpp"
#include "models/parser.hpp"
#include "serialization/formats.hpp"
#include "serialization/serializer_impl.hpp"

// clang-format off
namespace serialization {
using namespace lahuta;

template<>
struct Serializer<fmt::binary, ModelParserResult> {
  using Record = ModelParserResult;
  static constexpr std::size_t CoordStride = CoordinateStride;

  static std::size_t coords_reserved_bytes(const Record& r) {
    if (r.coords.empty()) return 0;
    return r.coords.size() * CoordStride;
  }

  static std::size_t serialized_size(const Record &r) {
    const std::size_t header_bytes = sizeof(ModelPayloadHeader);
    const std::size_t seq_len      = r.sequence.size();
    const std::size_t coords_bytes = coords_reserved_bytes(r);
    const std::size_t taxonomy_len = r.metadata.ncbi_taxonomy_id.size();
    const std::size_t organism_len = r.metadata.organism_scientific.size();
    const std::size_t plddt_bytes  = r.plddt_per_residue.size() * sizeof(pLDDTCategory);
    const std::size_t dssp_bytes   = r.dssp_per_residue.size() * sizeof(DSSPAssignment);
    return header_bytes + seq_len + coords_bytes + taxonomy_len + organism_len + plddt_bytes + dssp_bytes;
  }

  static void serialize_into_buffer(const Record &r, char *dest) {
    const auto seq_len      = static_cast<std::uint32_t>(r.sequence.size());
    const auto coords_count = static_cast<std::uint32_t>(r.coords.size());
    const auto coords_bytes = static_cast<std::uint32_t>(r.coords.size() * CoordStride);
    const auto taxonomy_len = static_cast<std::uint32_t>(r.metadata.ncbi_taxonomy_id.size());
    const auto organism_len = static_cast<std::uint32_t>(r.metadata.organism_scientific.size());
    const auto plddt_bytes  = static_cast<std::uint32_t>(r.plddt_per_residue.size() * sizeof(pLDDTCategory));
    const auto dssp_bytes   = static_cast<std::uint32_t>(r.dssp_per_residue.size()  * sizeof(DSSPAssignment));

    ModelPayloadHeader header{};
    header.magic           = ModelPayloadMagic;
    header.version         = ModelPayloadVersion;
    header.flags           = 0;
    header.sequence_length = seq_len;
    header.coord_count     = coords_count;

    std::uint32_t cursor   = static_cast<std::uint32_t>(sizeof(ModelPayloadHeader));
    const auto total_bytes = static_cast<std::uint32_t>(serialized_size(r));

    auto assign_field = [&](ModelPayloadField field, const char* data, std::uint32_t len) {
      auto& slice = header.slices[static_cast<std::size_t>(field)];
      slice.length = len;
      if (len == 0) {
        slice.offset = 0;
        return;
      }
      slice.offset = cursor;
      std::memcpy(dest + cursor, data, len);
      cursor += len;
    };

    assign_field(ModelPayloadField::Sequence, r.sequence.data(), seq_len);

    auto& coord_slice = header.slices[static_cast<std::size_t>(ModelPayloadField::Coordinates)];
    const auto reserved_coords_bytes = static_cast<std::uint32_t>(coords_reserved_bytes(r));
    if (reserved_coords_bytes == 0) {
      coord_slice.offset = 0;
      coord_slice.length = 0;
    } else {
      coord_slice.offset = cursor;
      coord_slice.length = coords_bytes;

      for (std::size_t i = 0; i < r.coords.size(); ++i) {
        float triple[3] = {static_cast<float>(r.coords[i].x), static_cast<float>(r.coords[i].y), static_cast<float>(r.coords[i].z)};
        std::memcpy(dest + coord_slice.offset + static_cast<std::uint32_t>(i * CoordStride), triple, sizeof(triple));
      }
      cursor += reserved_coords_bytes;
    }

    assign_field(ModelPayloadField::TaxonomyId,         r.metadata.ncbi_taxonomy_id.data(),    taxonomy_len);
    assign_field(ModelPayloadField::OrganismScientific, r.metadata.organism_scientific.data(), organism_len);
    assign_field(ModelPayloadField::PlddtCategories, reinterpret_cast<const char*>(r.plddt_per_residue.data()), plddt_bytes);
    assign_field(ModelPayloadField::DsspAssignments, reinterpret_cast<const char*>(r.dssp_per_residue.data()),  dssp_bytes);

    if (cursor != total_bytes) throw std::runtime_error("Model payload serialization size mismatch");

    std::memcpy(dest, &header, sizeof(header));
  }

  static std::string serialize(const Record &r) {
    std::string out(serialized_size(r), '\0');
    serialize_into_buffer(r, out.data());
    return out;
  }

  static Record deserialize(const char *buf, std::size_t n) {
    ModelPayloadView view(std::string_view(buf, n));
    Record r;

    auto seq = view.sequence();
    r.sequence.assign(seq.begin(), seq.end());

    auto taxonomy = view.taxonomy_id();
    r.metadata.ncbi_taxonomy_id.assign(taxonomy.begin(), taxonomy.end());

    auto organism = view.organism_scientific();
    r.metadata.organism_scientific.assign(organism.begin(), organism.end());

    auto coords_bytes = view.coords_bytes();
    const std::size_t coord_count = view.coord_count();
    if (coord_count > 0) {
      if (coords_bytes.size() != coord_count * CoordStride) throw std::runtime_error("Corrupted coordinate payload");
      r.coords.resize(coord_count);
      for (std::size_t i = 0; i < coord_count; ++i) {
        float triple[3];
        std::memcpy(triple, coords_bytes.data() + i * sizeof(triple), sizeof(triple));
        r.coords[i].x = static_cast<double>(triple[0]);
        r.coords[i].y = static_cast<double>(triple[1]);
        r.coords[i].z = static_cast<double>(triple[2]);
      }
    } else {
      r.coords.clear();
    }

    auto plddt = view.plddt_bytes();
    if (!plddt.empty()) {
      if (plddt.size() % sizeof(pLDDTCategory) != 0) throw std::runtime_error("Corrupted pLDDT payload");
      const std::size_t len = plddt.size() / sizeof(pLDDTCategory);
      r.plddt_per_residue.resize(len);
      std::memcpy(r.plddt_per_residue.data(), plddt.data(), plddt.size());
    } else {
      r.plddt_per_residue.clear();
    }

    auto dssp = view.dssp_bytes();
    if (!dssp.empty()) {
      if (dssp.size() % sizeof(DSSPAssignment) != 0) throw std::runtime_error("Corrupted DSSP payload");
      const std::size_t len = dssp.size() / sizeof(DSSPAssignment);
      r.dssp_per_residue.resize(len);
      std::memcpy(r.dssp_per_residue.data(), dssp.data(), dssp.size());
    } else {
      r.dssp_per_residue.clear();
    }

    return r;
  }

  static Record deserialize(const std::string &s) { return deserialize(s.data(), s.size()); }
};

// Binary serialization for analysis::system::ModelRecord
template<>
struct Serializer<fmt::binary, analysis::system::ModelRecord> {
  using Rec = analysis::system::ModelRecord;

  static std::size_t serialized_size(const Rec& r) {
    const uint32_t path_len  = static_cast<uint32_t>(r.file_path.size());
    const uint32_t blob_len  = static_cast<uint32_t>(Serializer<fmt::binary, ModelParserResult>::serialized_size(r.data));
    return 1 + sizeof(uint32_t) + path_len + sizeof(uint32_t) + blob_len;
  }

  // Serialize into an existing buffer, resizing it exactly once, no intermmediate allocations.
  static void serialize_into(const Rec& r, std::string& out) {
    out.resize(serialized_size(r));
    serialize_into_buffer(r, out.data(), out.size());
  }

  static void serialize_into_buffer(const Rec& r, char* dest, std::size_t buffer_size) {
    const auto total = serialized_size(r);
    if (buffer_size < total) throw std::runtime_error("ModelRecord buffer too small");

    char* p = dest;
    const uint32_t path_len  = static_cast<uint32_t>(r.file_path.size());
    const uint32_t blob_len  = static_cast<uint32_t>(Serializer<fmt::binary, ModelParserResult>::serialized_size(r.data));

    *p++ = static_cast<char>(r.success ? 1 : 0);
    std::memcpy(p, &path_len, sizeof(path_len));
    p += sizeof(path_len);
    if (path_len) {
      std::memcpy(p, r.file_path.data(), path_len);
      p += path_len;
    }
    std::memcpy(p, &blob_len, sizeof(blob_len));
    p += sizeof(blob_len);
    Serializer<fmt::binary, ModelParserResult>::serialize_into_buffer(r.data, p);
  }

  static std::string serialize(const Rec &r) {
    std::string out;
    serialize_into(r, out);
    return out;
  }

  // FIX: we do not perform any validation (e.g. check size, endianness, etc.)
  static Rec deserialize(const char *data, std::size_t size) {
    Rec r;
    size_t off = 0;
    r.success = static_cast<unsigned char>(data[off++]);
    uint32_t path_len;
    std::memcpy(&path_len, data + off, sizeof(path_len));
    off += sizeof(path_len);
    r.file_path.assign(data + off, path_len);
    off += path_len;
    // model blob
    uint32_t blob_len;
    std::memcpy(&blob_len, data + off, sizeof(blob_len));
    off += sizeof(blob_len);
    r.data = Serializer<fmt::binary, ModelParserResult>::deserialize(data + off, blob_len);
    return r;
  }
  static Rec deserialize(const std::string &s) { return deserialize(s.data(), s.size()); }
};

} // namespace serialization

#endif // LAHUTA_SERIALIZATION_MODEL_SERIALIZER_HPP 
