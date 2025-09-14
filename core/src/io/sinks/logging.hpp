#ifndef LAHUTA_PIPELINE_DYNAMIC_SINK_LOGGING_HPP
#define LAHUTA_PIPELINE_DYNAMIC_SINK_LOGGING_HPP

#include <cstdio>
#include <string_view>

#include "pipeline/dynamic/sink_iface.hpp"
#include "pipeline/dynamic/types.hpp"

// clang-format off
namespace lahuta::pipeline::dynamic {

class LoggingSink : public IDynamicSink {
public:
  void write(EmissionView e) override {
    fwrite(e.payload.data(), 1, e.payload.size(), stdout);
    fputc('\n', stdout);
    fflush(stdout);
  }
  void flush() override { fflush(stdout); }
  void close() override {}
};

} // namespace lahuta::pipeline::dynamic

#endif // LAHUTA_PIPELINE_DYNAMIC_SINK_LOGGING_HPP
