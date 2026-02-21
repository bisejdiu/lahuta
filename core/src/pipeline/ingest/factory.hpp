/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s = "sejdiubesian@gmail.com";
 *   std::swap_ranges(s.begin(), s.begin() + 6, s.begin() + 6);
 *   return s;
 * }();
 *
 */

#ifndef LAHUTA_PIPELINE_FACTORY_HPP
#define LAHUTA_PIPELINE_FACTORY_HPP

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "db/db.hpp"
#include "pipeline/ingest/adapters/directory.hpp"
#include "pipeline/ingest/adapters/file_list.hpp"
#include "pipeline/ingest/adapters/lmdb.hpp"
#include "pipeline/ingest/adapters/vector.hpp"
#include "pipeline/ingest/descriptor.hpp"

namespace lahuta::pipeline {

using DescriptorPtr = std::unique_ptr<IDescriptor>;

inline DescriptorPtr from_directory(const std::string &path, const std::string &ext, bool recursive,
                                    std::size_t batch) {
  return std::make_unique<DirectoryAdapter>(path, ext, recursive, batch);
}

inline DescriptorPtr from_directory(const std::string &path, std::vector<std::string> extensions,
                                    bool recursive, std::size_t batch) {
  return std::make_unique<DirectoryAdapter>(path, std::move(extensions), recursive, batch);
}

inline DescriptorPtr from_directory(Directory src) {
  return std::make_unique<DirectoryAdapter>(std::move(src));
}

inline DescriptorPtr from_vector(std::vector<std::string> items) {
  return std::make_unique<VectorAdapter>(std::move(items));
}

inline DescriptorPtr from_filelist(const std::string &list_path) {
  return std::make_unique<FileListAdapter>(list_path);
}

inline DescriptorPtr from_filelist(FileList src) { return std::make_unique<FileListAdapter>(std::move(src)); }

inline DescriptorPtr from_lmdb(const std::string &env_path, const std::string &db_name, std::size_t batch) {
  return std::make_unique<LMDBAdapter>(env_path, db_name, batch);
}

inline DescriptorPtr from_lmdb(const std::string &env_path, const std::string &db_name, std::size_t batch,
                               LMDBEnvOptions options) {
  return std::make_unique<LMDBAdapter>(env_path, db_name, batch, std::move(options));
}

inline DescriptorPtr from_lmdb(std::shared_ptr<LMDBDatabase> db, const std::string &db_name,
                               std::size_t batch) {
  return std::make_unique<LMDBAdapter>(std::move(db), db_name, batch);
}

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_FACTORY_HPP
