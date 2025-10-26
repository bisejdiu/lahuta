#ifndef LAHUTA_SERIALIZATION_MODEL_SERIALIZER_HPP
#define LAHUTA_SERIALIZATION_MODEL_SERIALIZER_HPP

#include <cstring>

#include "analysis/system/records.hpp"
#include "models/parser.hpp"
#include "serialization/formats.hpp"
#include "serialization/serializer_impl.hpp"

// clang-format off
namespace serialization {
using namespace lahuta;

template<>
struct Serializer<fmt::binary, ModelParserResult> {
  using Record = ModelParserResult;

  static std::size_t serialized_size(const Record &r) {
    const std::size_t seq_len       = r.sequence.size();
    const std::size_t coords_bytes  = r.coords.size() * 3 * sizeof(float);
    const std::size_t taxonomy_len  = r.metadata.ncbi_taxonomy_id.size();
    const std::size_t organism_len  = r.metadata.organism_scientific.size();
    const std::size_t plddt_bytes   = r.plddt_per_residue.size() * sizeof(pLDDTCategory);
    // Layout: [seq_len][n_points][sequence][coords][taxonomy_len][organism_len][taxonomy][organism][plddt_len][plddt_bytes]
    return sizeof(uint32_t) * 2 + seq_len + coords_bytes + sizeof(uint32_t) * 2 +
           taxonomy_len + organism_len + sizeof(uint32_t) + plddt_bytes;
  }

  static void serialize_into_buffer(const Record &r, char *dest) {
    const uint32_t seq_len      = static_cast<uint32_t>(r.sequence.size());
    const uint32_t n_points     = static_cast<uint32_t>(r.coords.size());
    const uint32_t taxonomy_len = static_cast<uint32_t>(r.metadata.ncbi_taxonomy_id.size());
    const uint32_t organism_len = static_cast<uint32_t>(r.metadata.organism_scientific.size());
    const uint32_t plddt_len    = static_cast<uint32_t>(r.plddt_per_residue.size());

    char *p = dest;
    std::memcpy(p, &seq_len, sizeof(seq_len));
    p += sizeof(seq_len);
    std::memcpy(p, &n_points, sizeof(n_points));
    p += sizeof(n_points);

    if (seq_len) {
      std::memcpy(p, r.sequence.data(), seq_len);
    }
    p += seq_len;

    for (const auto &v : r.coords) {
      float f[3] = {
          static_cast<float>(v.x),
          static_cast<float>(v.y),
          static_cast<float>(v.z)};
      std::memcpy(p, f, sizeof(f));
      p += sizeof(f);
    }

    std::memcpy(p, &taxonomy_len, sizeof(taxonomy_len));
    p += sizeof(taxonomy_len);
    std::memcpy(p, &organism_len, sizeof(organism_len));
    p += sizeof(organism_len);

    if (taxonomy_len) {
      std::memcpy(p, r.metadata.ncbi_taxonomy_id.data(), taxonomy_len);
    }
    p += taxonomy_len;

    if (organism_len) {
      std::memcpy(p, r.metadata.organism_scientific.data(), organism_len);
    }
    p += organism_len;

    std::memcpy(p, &plddt_len, sizeof(plddt_len));
    p += sizeof(plddt_len);
    if (plddt_len) {
      const size_t plddt_bytes = static_cast<size_t>(plddt_len) * sizeof(pLDDTCategory);
      std::memcpy(p, r.plddt_per_residue.data(), plddt_bytes);
      p += plddt_bytes;
    }
  }

  static std::string serialize(const Record &r) {
    std::string out(serialized_size(r), '\0');
    serialize_into_buffer(r, out.data());
    return out;
  }

  static Record deserialize(const char *buf, std::size_t n) {
    Record r;

    if (n < sizeof(uint32_t) * 2) {
      throw std::runtime_error("Corrupted data in deserialization");
    }

    const char *p = buf;
    uint32_t seq_len, n_points;
    std::memcpy(&seq_len, p, sizeof(seq_len));
    p += sizeof(seq_len);
    std::memcpy(&n_points, p, sizeof(n_points));
    p += sizeof(n_points);

    const std::size_t coords_bytes = static_cast<std::size_t>(n_points) * 3 * sizeof(float);
    const std::size_t base_size = sizeof(uint32_t) * 2 + seq_len + coords_bytes;
    if (n < base_size) {
      throw std::runtime_error("Corrupted data in deserialization");
    }

    r.sequence.assign(p, seq_len);
    p += seq_len;

    r.coords.resize(n_points);
    for (auto &v : r.coords) {
      float f[3];
      std::memcpy(f, p, sizeof(f));
      v.x = double(f[0]);
      v.y = double(f[1]);
      v.z = double(f[2]);
      p += sizeof(f);
    }

    if (n > static_cast<std::size_t>(p - buf)) {
      std::size_t remaining = n - static_cast<std::size_t>(p - buf);
      if (remaining < sizeof(uint32_t) * 2) {
        throw std::runtime_error("Corrupted data in deserialization (taxonomy fields)");
      }

      uint32_t taxonomy_len = 0;
      uint32_t organism_len = 0;
      std::memcpy(&taxonomy_len, p, sizeof(taxonomy_len));
      p += sizeof(taxonomy_len);
      std::memcpy(&organism_len, p, sizeof(organism_len));
      p += sizeof(organism_len);
      remaining -= sizeof(uint32_t) * 2;

      if (remaining < static_cast<std::size_t>(taxonomy_len) + static_cast<std::size_t>(organism_len)) {
        throw std::runtime_error("Corrupted data in deserialization (taxonomy payload)");
      }

      r.metadata.ncbi_taxonomy_id.assign(p, taxonomy_len);
      p += taxonomy_len;
      r.metadata.organism_scientific.assign(p, organism_len);
      p += organism_len;
      remaining -= static_cast<std::size_t>(taxonomy_len) + static_cast<std::size_t>(organism_len);

      if (remaining >= sizeof(uint32_t)) {
        uint32_t plddt_len = 0;
        std::memcpy(&plddt_len, p, sizeof(plddt_len));
        p += sizeof(plddt_len);
        remaining -= sizeof(plddt_len);
        const size_t plddt_bytes = static_cast<size_t>(plddt_len) * sizeof(pLDDTCategory);
        if (remaining < plddt_bytes) {
          throw std::runtime_error("Corrupted data in deserialization (pLDDT payload)");
        }
        r.plddt_per_residue.resize(plddt_len);
        if (plddt_bytes) {
          std::memcpy(r.plddt_per_residue.data(), p, plddt_bytes);
          p += plddt_bytes;
          remaining -= plddt_bytes;
        } else {
          r.plddt_per_residue.clear();
        }
      } else {
        r.plddt_per_residue.clear();
      }
    } else {
      r.metadata.ncbi_taxonomy_id.clear();
      r.metadata.organism_scientific.clear();
      r.plddt_per_residue.clear();
    }
    return r;
  }
  static Record deserialize(const std::string &s) { return deserialize(s.data(), s.size()); }
};

// Binary serialization for analysis::system::ModelRecord
template<>
struct Serializer<fmt::binary, analysis::system::ModelRecord> {
  using Rec = analysis::system::ModelRecord;

  // Serialize into an existing buffer, resizing it exactly once, no intermmediate allocations.
  static void serialize_into(const Rec& r, std::string& out) {
    const uint32_t path_len  = static_cast<uint32_t>(r.file_path.size());
    const uint32_t blob_len  = static_cast<uint32_t>(
        Serializer<fmt::binary, ModelParserResult>::serialized_size(r.data));

    const std::size_t total = 1 + sizeof(uint32_t) + path_len + sizeof(uint32_t) + blob_len;
    out.resize(total);

    char* p = out.data();
    // success byte
    *p++ = static_cast<char>(r.success ? 1 : 0);
    // path length
    std::memcpy(p, &path_len, sizeof(path_len));
    p += sizeof(path_len);
    // path
    if (path_len) {
      std::memcpy(p, r.file_path.data(), path_len);
      p += path_len;
    }
    // blob length
    std::memcpy(p, &blob_len, sizeof(blob_len));
    p += sizeof(blob_len);
    // ModelParserResult blob = [seq_len][n_points][sequence][coords(float3)*n][taxonomy_len][organism_len][taxonomy][organism]
    Serializer<fmt::binary, ModelParserResult>::serialize_into_buffer(r.data, p);
    p += blob_len;
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
