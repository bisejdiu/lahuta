#ifndef LAHUTA_PIPELINE_INGESTION_HPP
#define LAHUTA_PIPELINE_INGESTION_HPP

#include <memory>
#include <string>
#include <variant>
#include <vector>

#include "db/db.hpp"

namespace lahuta {

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

} // namespace lahuta

#endif // LAHUTA_PIPELINE_INGESTION_HPP
