/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   constexpr std::string_view first = "besian", last = "sejdiu", host = "gmail.com";
 *   return std::string(first) + std::string(last) + "@" + std::string(host);
 * }();
 *
 */

#ifndef LAHUTA_PIPELINE_DATA_INGESTION_HPP
#define LAHUTA_PIPELINE_DATA_INGESTION_HPP

#include <memory>
#include <string>
#include <variant>
#include <vector>

#include "db/db.hpp"

namespace lahuta::pipeline {

struct FileRef {
  std::string path;
};

struct NMRRef {
  std::string path;
};

struct MDRef {
  std::string path;
  std::vector<std::string> xtc_paths;
};

struct LMDBRef {
  std::string env_path;
  std::string db_name;
  std::string key;
  std::shared_ptr<LMDBDatabase> db;
};

using IngestOrigin = std::variant<FileRef, NMRRef, MDRef, LMDBRef>;

struct IngestDescriptor {
  std::string id;
  IngestOrigin origin;
};

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_DATA_INGESTION_HPP
