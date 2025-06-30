#ifndef LAHUTA_SERIALIZATION_MODEL_SERIALIZER_HPP
#define LAHUTA_SERIALIZATION_MODEL_SERIALIZER_HPP

#include "models/parser.hpp"
#include "serialization/serializer_impl.hpp"
#include "tasks/model_db_writer.hpp"
#include <cstring>
#include <serialization/formats.hpp>

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

template<>
struct Serializer<fmt::binary, tasks::ModelWriteTask::result_type> {
  using Rec = tasks::ModelWriteTask::result_type;
  static std::string serialize(const Rec &r) {
    std::string out;

    auto model_blob = Serializer<fmt::binary, ModelParserResult>::serialize(r.data);
    uint32_t path_len = static_cast<uint32_t>(r.file_path.size());
    uint32_t blob_len = static_cast<uint32_t>(model_blob.size());
    out.reserve(1 + sizeof(path_len) + path_len + sizeof(blob_len) + blob_len);

    out.push_back(r.success ? 1 : 0);
    out.append(reinterpret_cast<const char *>(&path_len), sizeof(path_len));
    out.append(r.file_path);
    out.append(reinterpret_cast<const char *>(&blob_len), sizeof(blob_len));
    out.append(model_blob);
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
