#ifndef LAHUTA_PIPELINE_DYNAMIC_SOURCES_HPP
#define LAHUTA_PIPELINE_DYNAMIC_SOURCES_HPP

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "db/db.hpp"
#include "sources/adapters/directory.hpp"
#include "sources/adapters/file_list.hpp"
#include "sources/adapters/lmdb.hpp"
#include "sources/adapters/vector.hpp"
#include "sources/descriptor.hpp"

// clang-format off
namespace lahuta::pipeline::dynamic::sources_factory {

using DescriptorPtr = std::unique_ptr<sources::IDescriptor>;

inline DescriptorPtr from_directory(const std::string& path, const std::string& ext, bool recursive, std::size_t batch) {
  return std::make_unique<sources::DirectoryAdapter>(path, ext, recursive, batch);
}

inline DescriptorPtr from_directory(sources::Directory src) {
  return std::make_unique<sources::DirectoryAdapter>(std::move(src));
}

inline DescriptorPtr from_vector(std::vector<std::string> items) {
  return std::make_unique<sources::VectorAdapter>(std::move(items));
}

inline DescriptorPtr from_filelist(const std::string& list_path) {
  return std::make_unique<sources::FileListAdapter>(list_path);
}

inline DescriptorPtr from_filelist(sources::FileList src) {
  return std::make_unique<sources::FileListAdapter>(std::move(src));
}

inline DescriptorPtr from_lmdb(const std::string& env_path, const std::string& db_name, std::size_t batch) {
  return std::make_unique<sources::LMDBAdapter>(env_path, db_name, batch);
}

inline DescriptorPtr from_lmdb(std::shared_ptr<LMDBDatabase> db, const std::string& db_name, std::size_t batch) {
  return std::make_unique<sources::LMDBAdapter>(std::move(db), db_name, batch);
}

} // namespace lahuta::pipeline::dynamic::sources_factory

#endif // LAHUTA_PIPELINE_DYNAMIC_SOURCES_HPP
