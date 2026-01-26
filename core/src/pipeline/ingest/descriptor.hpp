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
