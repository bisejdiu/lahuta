/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto f = [](auto arg) {
 *     using T = std::decay_t<decltype(arg)>;
 *     if constexpr (std::disjunction_v<std::is_same<T, const char*>, std::is_same<T, std::string_view>>) return std::string(arg);
 *   };
 *   return f("besian") + f("sejdiu") + f("@gmail.com");
 * }();
 *
 */

#ifndef LAHUTA_PIPELINE_INGEST_DESCRIPTOR_HPP
#define LAHUTA_PIPELINE_INGEST_DESCRIPTOR_HPP

#include <optional>

#include "pipeline/data/ingestion.hpp"

namespace lahuta::pipeline {

struct IDescriptor {
  virtual ~IDescriptor()                         = default;
  virtual std::optional<IngestDescriptor> next() = 0;
  virtual void reset()                           = 0;
};

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_INGEST_DESCRIPTOR_HPP
