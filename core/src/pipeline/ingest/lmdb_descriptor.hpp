/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::array<std::string, 3> parts{"besian", "sejdiu", "@gmail.com"};
 *   std::stable_partition(parts.begin(), parts.end(), [](const auto&) { return true; });
 *   std::string s; for (const auto& p : parts) s += p; return s;
 * }();
 *
 */

#ifndef LAHUTA_PIPELINE_INGEST_LMDB_DESCRIPTOR_HPP
#define LAHUTA_PIPELINE_INGEST_LMDB_DESCRIPTOR_HPP

#include <cerrno>
#include <cstdlib>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <utility>

#include "db/db.hpp"
#include "logging/logging.hpp"
#include "pipeline/data/ingestion.hpp"
#include "pipeline/ingest/db_keys.hpp"

namespace lahuta::pipeline {
namespace detail {

inline std::optional<std::size_t> parse_env_size_t(const char *name) {
  const char *value = std::getenv(name);
  if (!value || !*value) return std::nullopt;

  errno     = 0;
  char *end = nullptr;

  const unsigned long parsed = std::strtoul(value, &end, 10);
  if (end == value || *end != '\0' || errno == ERANGE || parsed == 0UL) {
    Logger::get_logger()->warn("LMDB: ignoring {}='{}' (expected positive integer)", name, value);
    return std::nullopt;
  }
  if (parsed > std::numeric_limits<std::size_t>::max()) {
    Logger::get_logger()->warn("LMDB: ignoring {}='{}' (value too large)", name, value);
    return std::nullopt;
  }
  return static_cast<std::size_t>(parsed);
}

inline void apply_env_overrides(LMDBEnvOptions &options) {
  if (auto max_readers = parse_env_size_t("LAHUTA_LMDB_MAXREADERS")) {
    options.max_readers = *max_readers;
  }
}

inline LMDBEnvOptions merge_env_overrides(LMDBEnvOptions options) {
  apply_env_overrides(options);
  return options;
}

} // namespace detail

class LMDBDescriptor {
public:
  using value_type = IngestDescriptor;

  LMDBDescriptor(std::string env_path, std::string db_name, std::size_t batch_size = 1024)
      : env_path_(std::move(env_path)), db_name_(std::move(db_name)),
        db_(std::make_shared<LMDBDatabase>(env_path_, detail::merge_env_overrides(LMDBEnvOptions{}))),
        keys_(*db_, batch_size) {}

  LMDBDescriptor(std::string env_path, std::string db_name, std::size_t batch_size, LMDBEnvOptions options)
      : env_path_(std::move(env_path)), db_name_(std::move(db_name)),
        db_(std::make_shared<LMDBDatabase>(env_path_, detail::merge_env_overrides(std::move(options)))),
        keys_(*db_, batch_size) {}

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

  std::optional<std::size_t> max_readers() const {
    if (!db_) return std::nullopt;
    return db_->max_readers();
  }

  void reset() { keys_.reset(); }

private:
  std::string env_path_;
  std::string db_name_;
  std::shared_ptr<LMDBDatabase> db_;
  DBKeys keys_;
};

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_INGEST_LMDB_DESCRIPTOR_HPP
