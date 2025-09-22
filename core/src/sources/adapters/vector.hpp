#ifndef LAHUTA_SOURCES_ADAPTERS_VECTOR_HPP
#define LAHUTA_SOURCES_ADAPTERS_VECTOR_HPP

#include <optional>
#include <string>
#include <vector>

#include "sources/descriptor.hpp"
#include "sources/in_memory.hpp"

// clang-format off
namespace lahuta::sources {

class VectorAdapter final : public IDescriptor {
public:
  explicit VectorAdapter(std::vector<std::string> items) : inner_(std::move(items)) {}
  explicit VectorAdapter(sources::InMemory<std::string> src) : inner_(std::move(src)) {}

  std::optional<IngestDescriptor> next() override {
    auto val_ptr = inner_.next();
    if (!val_ptr) return std::nullopt;
    IngestDescriptor desc;
    const auto &val = *(*val_ptr);
    desc.id = val;
    desc.origin = FileRef{val};
    return desc;
  }
  void reset() override { inner_.reset(); }

private:
  sources::InMemory<std::string> inner_;
};

} // namespace lahuta::sources

#endif // LAHUTA_SOURCES_ADAPTERS_VECTOR_HPP
