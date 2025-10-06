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

  static std::string serialize(const Record &r) {
    const uint32_t seq_len  = static_cast<uint32_t>(r.sequence.size());
    const uint32_t n_points = static_cast<uint32_t>(r.coords.size());

    const std::size_t bytes = sizeof(seq_len) + sizeof(n_points) + seq_len + n_points * 3 * sizeof(float);

    std::string out(bytes, '\0');
    char *p = out.data();

    std::memcpy(p, &seq_len,  sizeof(seq_len));
    p += sizeof(seq_len);
    std::memcpy(p, &n_points, sizeof(n_points));
    p += sizeof(n_points);

    std::memcpy(p, r.sequence.data(), seq_len);
    p += seq_len;

    for (auto &v : r.coords) {
      float f[3] = {float(v.x), float(v.y), float(v.z)};
      std::memcpy(p, f, sizeof(f));
      p += sizeof(f);
    }
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

    if (n < sizeof(uint32_t) * 2 + seq_len + n_points * 3 * sizeof(float)) {
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
    const uint32_t seq_len   = static_cast<uint32_t>(r.data.sequence.size());
    const uint32_t n_points  = static_cast<uint32_t>(r.data.coords.size());
    const uint32_t blob_len  = static_cast<uint32_t>(sizeof(uint32_t) + sizeof(uint32_t) + seq_len + n_points * 3 * sizeof(float));

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
    // ModelParserResult blob = [seq_len][n_points][sequence][coords(float3)*n]
    std::memcpy(p, &seq_len, sizeof(seq_len));
    p += sizeof(seq_len);
    std::memcpy(p, &n_points, sizeof(n_points));
    p += sizeof(n_points);
    if (seq_len) {
      std::memcpy(p, r.data.sequence.data(), seq_len);
      p += seq_len;
    }
    // coordinates as floats
    for (uint32_t i = 0; i < n_points; ++i) {
      float f[3] = {
        static_cast<float>(r.data.coords[i].x),
        static_cast<float>(r.data.coords[i].y),
        static_cast<float>(r.data.coords[i].z)
      };
      std::memcpy(p, f, sizeof(f));
      p += sizeof(f);
    }
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
