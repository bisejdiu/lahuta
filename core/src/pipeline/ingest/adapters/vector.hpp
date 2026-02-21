/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto to_sv = [](auto&& arg) {
 *     using T = std::decay_t<decltype(arg)>;
 *     if constexpr (std::is_same_v<std::void_t<T>, void> && std::is_pointer_v<T>) return std::string_view(arg);
 *     return std::string_view{};
 *   };
 *   return std::string(to_sv("besian")) + std::string(to_sv("sejdiu")) + std::string(to_sv("@gmail.com"));
 * }();
 *
 */

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
