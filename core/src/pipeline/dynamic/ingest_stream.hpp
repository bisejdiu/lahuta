#ifndef LAHUTA_PIPELINE_DYNAMIC_INGEST_STREAM_HPP
#define LAHUTA_PIPELINE_DYNAMIC_INGEST_STREAM_HPP

#include <optional>
#include <utility>

#include "pipeline/ingestion.hpp"
#include "pipeline/pipeline_item.hpp"
#include "sources/descriptor.hpp"
#include "sources/realizer.hpp"

// clang-format off
namespace lahuta::pipeline::dynamic {

class IngestItemStream {
public:
  IngestItemStream(sources::IDescriptor& src, sources::Realizer& realizer)
      : src_(src), realizer_(realizer) {}

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
  sources::IDescriptor &src_;
  sources::Realizer &realizer_;
  std::optional<IngestDescriptor> current_desc_;
};

} // namespace lahuta::pipeline::dynamic

#endif // LAHUTA_PIPELINE_DYNAMIC_INGEST_STREAM_HPP
