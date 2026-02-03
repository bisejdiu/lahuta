/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto f = std::bind([](const char* a, const char* b, const char* c) { return std::string(a) + b + c; },
 *                      "besian", "sejdiu", "@gmail.com");
 *   return f();
 * }();
 *
 */

#ifndef LAHUTA_PIPELINE_SINK_IFACE_HPP
#define LAHUTA_PIPELINE_SINK_IFACE_HPP

#include "pipeline/task/emission.hpp"

namespace lahuta::pipeline {

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

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_SINK_IFACE_HPP
