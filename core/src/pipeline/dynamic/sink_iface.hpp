#ifndef LAHUTA_PIPELINE_DYNAMIC_SINK_IFACE_HPP
#define LAHUTA_PIPELINE_DYNAMIC_SINK_IFACE_HPP

#include "pipeline/dynamic/types.hpp"

// clang-format off
namespace lahuta::pipeline::dynamic {

class IDynamicSink {
public:
  virtual ~IDynamicSink() = default;

  // Writer threads (potentially more than one per sink) call write(e).
  // Implementations must be thread-safe and must not retain e.payload beyond this call.
  virtual void write(EmissionView e) = 0;

  // Optional hooks
  virtual void flush() {}
  virtual void close() {}
};

} // namespace lahuta::pipeline::dynamic

#endif // LAHUTA_PIPELINE_DYNAMIC_SINK_IFACE_HPP
