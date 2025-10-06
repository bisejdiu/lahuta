#ifndef LAHUTA_SOURCES_ADAPTERS_FILE_LIST_HPP
#define LAHUTA_SOURCES_ADAPTERS_FILE_LIST_HPP

#include <optional>
#include <string>

#include "sources/descriptor.hpp"
#include "sources/file_list.hpp"

// clang-format off
namespace lahuta::sources {

class FileListAdapter final : public IDescriptor {
public:
  explicit FileListAdapter(const std::string &list_path) : inner_(list_path) {}
  explicit FileListAdapter(sources::FileList src) : inner_(std::move(src)) {}

  std::optional<IngestDescriptor> next() override {
    auto path = inner_.next();
    if (!path) return std::nullopt;
    IngestDescriptor desc;
    desc.id = *path;
    desc.origin = FileRef{*path};
    return desc;
  }
  void reset() override { inner_.reset(); }

private:
  sources::FileList inner_;
};

} // namespace lahuta::sources

#endif // LAHUTA_SOURCES_ADAPTERS_FILE_LIST_HPP
