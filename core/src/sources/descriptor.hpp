#ifndef LAHUTA_SOURCES_DESCRIPTOR_HPP
#define LAHUTA_SOURCES_DESCRIPTOR_HPP

#include <optional>

#include "pipeline/ingestion.hpp"

// clang-format off
namespace lahuta::sources {

struct IDescriptor {
  virtual ~IDescriptor() = default;
  virtual std::optional<IngestDescriptor> next() = 0;
  virtual void reset() = 0;
};

} // namespace lahuta::sources

#endif // LAHUTA_SOURCES_DESCRIPTOR_HPP
