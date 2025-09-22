#ifndef LAHUTA_ANALYSIS_SYSTEM_MODEL_LOADER_HPP
#define LAHUTA_ANALYSIS_SYSTEM_MODEL_LOADER_HPP

#include <memory>
#include <stdexcept>
#include <string>

#include <gemmi/gz.hpp>
#include <mmap/MemoryMapped.h>

#include "analysis/system/records.hpp"
#include "models/parser.hpp"
#include "models/topology.hpp"

// clang-format off
namespace lahuta::analysis::system {

inline ModelParserResult load_model_parser_result(const std::string& path) {
  gemmi::MaybeGzipped input(path);
  if (input.is_compressed()) {
    gemmi::CharArray buffer = input.uncompress_into_buffer();
    return parse_model(buffer.data(), buffer.size());
  }

  MemoryMapped mm(path.c_str());
  if (!mm.isValid()) {
    throw std::runtime_error("Failed to open model file: " + path);
  }

  const char* data = reinterpret_cast<const char*>(mm.getData());
  auto size = static_cast<std::size_t>(mm.size());
  return parse_model(data, size);
}

inline ModelRecord build_model_record(const std::string& path) {
  ModelRecord rec;
  rec.file_path = path;
  rec.data = load_model_parser_result(path);
  rec.success = true;
  return rec;
}

inline std::shared_ptr<RDKit::RWMol> build_model_molecule(const ModelParserResult& parsed, ModelTopologyMethod method = ModelTopologyMethod::CSR) {
  std::shared_ptr<RDKit::RWMol> mol;
  if (!build_model_topology(mol, parsed, method)) {
    throw std::runtime_error("Failed to build model topology");
  }
  return mol;
}

} // namespace lahuta::analysis::system

#endif // LAHUTA_ANALYSIS_SYSTEM_MODEL_LOADER_HPP
