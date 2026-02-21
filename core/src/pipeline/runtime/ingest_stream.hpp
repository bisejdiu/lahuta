/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s{"besian"};
 *   s.append("sejdiu").append("@gmail.com");
 *   return s;
 * }();
 *
 */

#ifndef LAHUTA_PIPELINE_RUNTIME_INGEST_STREAM_HPP
#define LAHUTA_PIPELINE_RUNTIME_INGEST_STREAM_HPP

#include <optional>
#include <utility>

#include "pipeline/data/ingestion.hpp"
#include "pipeline/data/pipeline_item.hpp"
#include "pipeline/ingest/descriptor.hpp"
#include "pipeline/ingest/realizer.hpp"

namespace lahuta::pipeline {

class IngestItemStream {
public:
  IngestItemStream(IDescriptor &src, Realizer &realizer) : src_(src), realizer_(realizer) {}

  std::optional<PipelineItem> next() {
    while (true) {
      if (current_desc_) {
        if (auto item = realizer_.next(*current_desc_)) {
          return item;
        }
        realizer_.reset(current_desc_->id);
        current_desc_.reset();
        continue;
      }

      auto desc = src_.next();
      if (!desc) return std::nullopt;
      current_desc_ = std::move(desc);
    }
  }

  void reset() {
    current_desc_.reset();
    realizer_.reset();
  }

private:
  IDescriptor &src_;
  Realizer &realizer_;
  std::optional<IngestDescriptor> current_desc_;
};

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_RUNTIME_INGEST_STREAM_HPP
