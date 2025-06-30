#ifndef LAHUTA_SERIALIZATION_TOPOLOGY_SERIALIZER_HPP
#define LAHUTA_SERIALIZATION_TOPOLOGY_SERIALIZER_HPP

#include "serialization/json.hpp"
#include "serialization/serializer_impl.hpp"
#include "tasks/topology_task.hpp"
#include <serialization/formats.hpp>
#include <sstream>

// clang-format off
namespace serialization {
using namespace lahuta;

using TopoRes = tasks::TopologyTask::result_type;

template<>
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

template<>
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

template<>
struct Serializer<fmt::binary, TopoRes> {
  using Record = TopoRes;
  static std::string serialize(const TopoRes &v) {

    std::string buffer;
    buffer.reserve(
      sizeof(size_t) + v.file_path.size() +
      sizeof(v.data.num_atoms) + sizeof(v.data.num_bonds) + sizeof(v.success));

    size_t path_len = v.file_path.size();
    buffer.append(reinterpret_cast<const char *>(&path_len), sizeof(path_len));
    buffer.append(v.file_path);

    buffer.append(reinterpret_cast<const char *>(&v.data.num_atoms), sizeof(v.data.num_atoms));
    buffer.append(reinterpret_cast<const char *>(&v.data.num_bonds), sizeof(v.data.num_bonds));
    buffer.append(reinterpret_cast<const char *>(&v.success), sizeof(v.success));

    return buffer;
  }

  static TopoRes deserialize(const char *data, std::size_t size) {
    TopoRes result;
    std::size_t offset = 0;

    std::size_t path_len;
    std::memcpy(&path_len, data + offset, sizeof(path_len));
    offset += sizeof(path_len);

    result.file_path.assign(data + offset, path_len);
    offset += path_len;

    std::memcpy(&result.data.num_atoms, data + offset, sizeof(result.data.num_atoms));
    offset += sizeof(result.data.num_atoms);
    std::memcpy(&result.data.num_bonds, data + offset, sizeof(result.data.num_bonds));
    offset += sizeof(result.data.num_bonds);
    std::memcpy(&result.success, data + offset, sizeof(result.success));
    return result;
  }

  static TopoRes deserialize(const std::string &buf) { return deserialize(buf.data(), buf.size()); }
};

} // namespace serialization

#endif // LAHUTA_SERIALIZATION_TOPOLOGY_SERIALIZER_HPP 
