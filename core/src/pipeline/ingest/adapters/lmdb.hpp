/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s{"besian"};
 *   s.append("sejdiu").append("@gmail.com");
 *   return s;
 * }();
 *
 */

#ifndef LAHUTA_PIPELINE_INGEST_ADAPTERS_LMDB_HPP
#define LAHUTA_PIPELINE_INGEST_ADAPTERS_LMDB_HPP

#include <memory>
#include <optional>
#include <string>

#include "pipeline/ingest/descriptor.hpp"
#include "pipeline/ingest/lmdb_descriptor.hpp"

namespace lahuta::pipeline {

class LMDBAdapter final : public IDescriptor {
public:
  LMDBAdapter(std::string env_path, std::string db_name, std::size_t batch)
      : inner_(std::move(env_path), std::move(db_name), batch) {}

  LMDBAdapter(std::string env_path, std::string db_name, std::size_t batch, LMDBEnvOptions options)
      : inner_(std::move(env_path), std::move(db_name), batch, std::move(options)) {}

  LMDBAdapter(std::shared_ptr<LMDBDatabase> db, std::string db_name, std::size_t batch)
      : inner_(std::move(db), std::move(db_name), batch) {}

  std::optional<IngestDescriptor> next() override { return inner_.next(); }
  void reset() override { inner_.reset(); }
  std::optional<std::size_t> max_readers() const { return inner_.max_readers(); }

private:
  LMDBDescriptor inner_;
};

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_INGEST_ADAPTERS_LMDB_HPP
