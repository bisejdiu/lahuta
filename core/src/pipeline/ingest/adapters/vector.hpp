#ifndef LAHUTA_PIPELINE_INGEST_ADAPTERS_VECTOR_HPP
#define LAHUTA_PIPELINE_INGEST_ADAPTERS_VECTOR_HPP

#include <optional>
#include <string>
#include <vector>

#include "pipeline/ingest/descriptor.hpp"
#include "pipeline/ingest/in_memory.hpp"

namespace lahuta::pipeline {

class VectorAdapter final : public IDescriptor {
public:
  explicit VectorAdapter(std::vector<std::string> items) : inner_(std::move(items)) {}
  explicit VectorAdapter(InMemory<std::string> src) : inner_(std::move(src)) {}

  std::optional<IngestDescriptor> next() override {
    auto val_ptr = inner_.next();
    if (!val_ptr) return std::nullopt;
    IngestDescriptor desc;
    const auto &val = *(*val_ptr);
    desc.id         = val;
    desc.origin     = FileRef{val};
    return desc;
  }
  void reset() override { inner_.reset(); }

private:
  InMemory<std::string> inner_;
};

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_INGEST_ADAPTERS_VECTOR_HPP
