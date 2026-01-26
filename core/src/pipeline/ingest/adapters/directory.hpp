#ifndef LAHUTA_PIPELINE_INGEST_ADAPTERS_DIRECTORY_HPP
#define LAHUTA_PIPELINE_INGEST_ADAPTERS_DIRECTORY_HPP

#include <initializer_list>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "pipeline/ingest/descriptor.hpp"
#include "pipeline/ingest/directory.hpp"

namespace lahuta::pipeline {

class DirectoryAdapter final : public IDescriptor {
public:
  explicit DirectoryAdapter(Directory src) : inner_(std::move(src)) {}
  DirectoryAdapter(const std::string &path, const std::string &ext, bool recursive, std::size_t batch)
      : inner_(path, ext, recursive, batch) {}
  DirectoryAdapter(const std::string &path, std::vector<std::string> extensions, bool recursive,
                   std::size_t batch)
      : inner_(path, std::move(extensions), recursive, batch) {}
  DirectoryAdapter(const std::string &path, std::initializer_list<std::string_view> extensions,
                   bool recursive, std::size_t batch)
      : inner_(path, extensions, recursive, batch) {}

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
  Directory inner_;
};

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_INGEST_ADAPTERS_DIRECTORY_HPP
