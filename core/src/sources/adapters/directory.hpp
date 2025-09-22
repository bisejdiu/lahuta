#ifndef LAHUTA_SOURCES_ADAPTERS_DIRECTORY_HPP
#define LAHUTA_SOURCES_ADAPTERS_DIRECTORY_HPP

#include <optional>
#include <string>

#include "sources/descriptor.hpp"
#include "sources/directory.hpp"

// clang-format off
namespace lahuta::sources {

class DirectoryAdapter final : public IDescriptor {
public:
  explicit DirectoryAdapter(sources::Directory src) : inner_(std::move(src)) {}
  DirectoryAdapter(const std::string &path, const std::string &ext, bool recursive, std::size_t batch)
      : inner_(path, ext, recursive, batch) {}

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
  sources::Directory inner_;
};

} // namespace lahuta::sources

#endif // LAHUTA_SOURCES_ADAPTERS_DIRECTORY_HPP
