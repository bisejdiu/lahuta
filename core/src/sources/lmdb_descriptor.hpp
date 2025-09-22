#ifndef LAHUTA_SOURCES_LMDB_DESCRIPTOR_HPP
#define LAHUTA_SOURCES_LMDB_DESCRIPTOR_HPP

#include <memory>
#include <optional>
#include <string>
#include <utility>

#include "pipeline/ingestion.hpp"
#include "sources/db_keys.hpp"

namespace lahuta::sources {

class LMDBDescriptor {
public:
  using value_type = IngestDescriptor;

  LMDBDescriptor(std::string env_path, std::string db_name, std::size_t batch_size = 1024)
      : env_path_(std::move(env_path)), db_name_(std::move(db_name)),
        db_(std::make_shared<LMDBDatabase>(env_path_)), keys_(*db_, batch_size) {}

  LMDBDescriptor(std::shared_ptr<LMDBDatabase> db, std::string db_name, std::size_t batch_size = 1024)
      : env_path_(), db_name_(std::move(db_name)), db_(std::move(db)), keys_(*db_, batch_size) {}

  std::optional<IngestDescriptor> next() {
    auto key = keys_.next();
    if (!key) return std::nullopt;

    IngestDescriptor desc;
    desc.id = env_path_.empty() ? *key : env_path_ + "::" + *key;
    LMDBRef ref;
    ref.env_path = env_path_;
    ref.db_name  = db_name_;
    ref.key      = *key;
    ref.db       = db_;
    desc.origin  = std::move(ref);
    return desc;
  }

  void reset() { keys_.reset(); }

private:
  std::string env_path_;
  std::string db_name_;
  std::shared_ptr<LMDBDatabase> db_;
  DBKeys keys_;
};

} // namespace lahuta::sources

#endif // LAHUTA_SOURCES_LMDB_DESCRIPTOR_HPP
