#ifndef LAHUTA_SERIALIZATION_SERIALIZER_HPP
#define LAHUTA_SERIALIZATION_SERIALIZER_HPP

#include "formats.hpp"
#include "json.hpp"
#include "models/parser.hpp"
#include "tasks/model_db_writer.hpp"
#include "tasks/topology_task.hpp"
#include <string>

//
// We are cramming everything into one file. Not good.
// Our TopologyTask::result_type contains only test data, so we can hold off for now,
// until we have an actual production task. - Besian, May 16, 2025
//

// clang-format off
namespace serialization {
using namespace lahuta;

// hard error when we forget to specialise
template <class FormatTag, class T, class = void> struct Serializer {
  static_assert(sizeof(FormatTag) == 0, "No Serializer<FormatTag, T> available. Define one.");
};

// forward const‑qualified T to the existing non‑const serializer
template <class FormatTag, class T>
struct Serializer<FormatTag, const T> : Serializer<FormatTag, T> {};

//////////////////////////////////////////////////////////////////////////////////////////////
using TopoRes = tasks::TopologyTask::result_type;
template <>
struct Serializer<fmt::json, TopoRes> {
  using Record = TopoRes;
  static std::string serialize(const TopoRes &v) {
    return JsonBuilder{}
        .key("file_path") .value(v.file_path)
        .key("num_atoms") .value(v.data.num_atoms)
        .key("num_bonds") .value(v.data.num_bonds)
        .key("success")   .value(v.success)
        .str();
  }

  static TopoRes deserialize(const std::string &s) {
    TopoRes out;
    JsonReader r{s};
    out.file_path       = r.get<std::string>("file_path");
    out.data.num_atoms  = r.get<size_t>     ("num_atoms");
    out.data.num_bonds  = r.get<size_t>     ("num_bonds");
    out.success         = r.get<bool>       ("success");
    return out;
  }
};

template <>
struct Serializer<fmt::text, TopoRes> {
  using Record = TopoRes;
  static std::string serialize(const TopoRes &v) {
    std::ostringstream oss;
    oss << (v.success ? "1" : "0") << " "
        << v.file_path << " "
        << v.data.num_atoms << " "
        << v.data.num_bonds;
    return oss.str();
  }

  static TopoRes deserialize(const std::string &data) {
    TopoRes result;
    std::istringstream iss(data);

    int success;
    iss >> success >> result.file_path >> result.data.num_atoms >> result.data.num_bonds;
    result.success = (success != 0);

    return result;
  }
};

template <>
struct Serializer<fmt::binary, TopoRes> {
  using Record = TopoRes;
  static std::string serialize(const TopoRes &v) {
    std::string buffer;

    buffer.reserve(v.file_path.size() + sizeof(int) * 2 + sizeof(bool) + sizeof(size_t));

    // file_path length and content
    size_t path_len = v.file_path.size();
    buffer.append(reinterpret_cast<const char *>(&path_len), sizeof(path_len));
    buffer.append(v.file_path);

    // other fields
    buffer.append(reinterpret_cast<const char *>(&v.data.num_atoms), sizeof(v.data.num_atoms));
    buffer.append(reinterpret_cast<const char *>(&v.data.num_bonds), sizeof(v.data.num_bonds));
    buffer.append(reinterpret_cast<const char *>(&v.success), sizeof(v.success));

    return buffer;
  }

  static TopoRes deserialize(const char *data, std::size_t size) {
    TopoRes result;
    std::size_t offset = 0;

    // file_path length
    std::size_t path_len;
    std::memcpy(&path_len, data + offset, sizeof(path_len));
    offset += sizeof(path_len);

    // file_path content
    result.file_path.assign(data + offset, path_len);
    offset += path_len;

    // other fields
    std::memcpy(&result.data.num_atoms, data + offset, sizeof(result.data.num_atoms));
    offset += sizeof(result.data.num_atoms);
    std::memcpy(&result.data.num_bonds, data + offset, sizeof(result.data.num_bonds));
    offset += sizeof(result.data.num_bonds);
    std::memcpy(&result.success, data + offset, sizeof(result.success));
    return result;
  }

  static TopoRes deserialize(const std::string &buf) { return deserialize(buf.data(), buf.size()); }
};

template <>
struct Serializer<fmt::binary, ModelParserResult> {
  using Record = ModelParserResult;

  static std::string serialize(const Record &r) {
    const uint32_t seq_len = static_cast<uint32_t>(r.sequence.size());
    const uint32_t n_points = static_cast<uint32_t>(r.coords.size());

    const std::size_t bytes = sizeof(seq_len) + sizeof(n_points) + seq_len + n_points * 3 * sizeof(float);

    std::string out(bytes, '\0');
    char *p = out.data();

    std::memcpy(p, &seq_len, sizeof(seq_len));
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

template <>
struct Serializer<fmt::binary, tasks::ModelWriteTask::result_type> {
  using Rec = tasks::ModelWriteTask::result_type;
  static std::string serialize(const Rec &r) {
    std::string out;

    auto model_blob = Serializer<fmt::binary, ModelParserResult>::serialize(r.data);
    uint32_t path_len = static_cast<uint32_t>(r.file_path.size());
    uint32_t blob_len = static_cast<uint32_t>(model_blob.size());
    out.reserve(1 + sizeof(path_len) + path_len + sizeof(blob_len) + blob_len);

    // success flag
    out.push_back(r.success ? 1 : 0);
    // file_path len + content
    out.append(reinterpret_cast<const char *>(&path_len), sizeof(path_len));
    out.append(r.file_path);
    // model blob len + content
    out.append(reinterpret_cast<const char *>(&blob_len), sizeof(blob_len));
    out.append(model_blob);
    return out;
  }
  static Rec deserialize(const char *data, std::size_t size) {
    Rec r;
    size_t off = 0;
    // success
    r.success = static_cast<unsigned char>(data[off++]);
    // file_path
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

#endif // LAHUTA_SERIALIZATION_SERIALIZER_HPP
