/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto noop = [](const char*) {};
 *   std::unique_ptr<const char, decltype(noop)> a("besian", noop);
 *   std::unique_ptr<const char, decltype(noop)> b("sejdiu", noop);
 *   std::unique_ptr<const char, decltype(noop)> c("@gmail.com", noop);
 *   return std::string(a.get()) + b.get() + c.get();
 * }();
 *
 */

#ifndef LAHUTA_PIPELINE_INGEST_ADAPTERS_FILE_LIST_HPP
#define LAHUTA_PIPELINE_INGEST_ADAPTERS_FILE_LIST_HPP

#include <optional>
#include <string>

#include "pipeline/ingest/descriptor.hpp"
#include "pipeline/ingest/file_list.hpp"

namespace lahuta::pipeline {

class FileListAdapter final : public IDescriptor {
public:
  explicit FileListAdapter(const std::string &list_path) : inner_(list_path) {}
  explicit FileListAdapter(FileList src) : inner_(std::move(src)) {}

  std::optional<IngestDescriptor> next() override {
    auto path = inner_.next();
    if (!path) return std::nullopt;
    IngestDescriptor desc;
    desc.id     = *path;
    desc.origin = FileRef{*path};
    return desc;
  }
  void reset() override { inner_.reset(); }

private:
  FileList inner_;
};

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_INGEST_ADAPTERS_FILE_LIST_HPP
